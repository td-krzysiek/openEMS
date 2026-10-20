// SPDX-License-Identifier: GPL-3.0-or-later
#include "mpi_transport.h"
#include <array>
#include <chrono>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>

namespace openems { namespace distributed {
namespace {
void Check(int error)
{
	if (error==MPI_SUCCESS) return;
	char message[MPI_MAX_ERROR_STRING]; int count=0;
	MPI_Error_string(error,message,&count);
	throw std::runtime_error(std::string("MPI domain runtime: ")+std::string(message,count));
}
void MainThread()
{
	int initialized=0,finalized=0,main_thread=0;
	MPI_Initialized(&initialized); MPI_Finalized(&finalized);
	if (!initialized || finalized) throw std::logic_error("MPI runtime is not active");
	Check(MPI_Is_thread_main(&main_thread));
	if (!main_thread) throw std::logic_error("MPI adapter requires the initializing thread");
}
struct Comm
{
	MPI_Comm value=MPI_COMM_NULL;
	explicit Comm(MPI_Comm source)
	{
		MainThread(); Check(MPI_Comm_dup(source,&value));
		Check(MPI_Comm_set_errhandler(value,MPI_ERRORS_RETURN));
	}
	~Comm()
	{
		int finalized=0; MPI_Finalized(&finalized);
		if (!finalized && value!=MPI_COMM_NULL) MPI_Comm_free(&value);
	}
	void Abort(const std::string& reason)
	{
		std::cerr<<"Distributed run failed: "<<reason<<std::endl;
		MPI_Abort(value,71);
		std::terminate();
	}
};

// MPI datatypes convert scalar representation; never transmit native structs.
typedef std::array<uint64_t,9> Header;
Header MakeHeader(const MessageId& id,size_t count)
{return {{0x4f454d53484c0001ULL,id.epoch.plan,id.epoch.step,uint64_t(id.epoch.phase),
	id.channel,id.source,id.destination,uint64_t(count),32}};}

class MPICompletion : public Completion
{
public:
	std::shared_ptr<Comm> comm;
	Header header,expected;
	HostBuffer buffer;
	std::array<MPI_Request,2> requests={{MPI_REQUEST_NULL,MPI_REQUEST_NULL}};
	bool receive=false,complete=false;
	explicit MPICompletion(std::shared_ptr<Comm> c):comm(c) {}
	~MPICompletion()
	{
		int finalized=0; MPI_Finalized(&finalized);
		if (!complete && !finalized && (requests[0]!=MPI_REQUEST_NULL || requests[1]!=MPI_REQUEST_NULL))
			comm->Abort("incomplete MPI buffer/request destroyed");
	}
	bool Ready() override
	{
		if (complete) return true;
		MainThread();
		std::array<MPI_Status,2> statuses; int done=0;
		Check(MPI_Testall(2,requests.data(),&done,statuses.data()));
		if (!done) return false;
		if (receive)
		{
			int header_count=0,payload_count=0;
			Check(MPI_Get_count(&statuses[0],MPI_UINT64_T,&header_count));
			Check(MPI_Get_count(&statuses[1],MPI_FLOAT,&payload_count));
			if (header_count!=9 || payload_count!=int(buffer->size()) || header!=expected)
				comm->Abort("wrong plan, epoch, phase, endpoint, scalar format or payload count");
		}
		complete=true; buffer.reset();
		return true;
	}
	void Wait() override
	{
		const auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(60);
		while (!Ready())
		{
			if (std::chrono::steady_clock::now()>deadline) comm->Abort("MPI exchange timeout");
			std::this_thread::yield();
		}
	}
};
}

struct MPIHaloTransport::Impl
{
	std::shared_ptr<Comm> comm;
	PartitionPlan plan;
	std::vector<uint32_t> ranks;
	int rank=0,size=0,tag_limit=0;
	std::vector<std::weak_ptr<MPICompletion>> pending;
	Impl(MPI_Comm source,const PartitionPlan& p,const std::vector<uint32_t>& r)
		:comm(std::make_shared<Comm>(source)),plan(p),ranks(r)
	{
		static_assert(sizeof(float)==4 && std::numeric_limits<float>::is_iec559,"FP32 wire contract required");
		Check(MPI_Comm_rank(comm->value,&rank)); Check(MPI_Comm_size(comm->value,&size));
		int* limit=nullptr,present=0;
		Check(MPI_Comm_get_attr(comm->value,MPI_TAG_UB,&limit,&present));
		bool valid=present && limit && ranks.size()==plan.Domains().size();
		if (valid)
		{
			tag_limit=*limit;
			for (auto target:ranks) valid=valid && target<uint32_t(size);
			for (const auto& recipe:plan.Halos())
				valid=valid && uint64_t(recipe.channel)*2+1<=uint64_t(tag_limit) && recipe.Samples()<=INT_MAX;
		}
		int local=valid?1:0,all=0;
		Check(MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,comm->value));
		if (!all) comm->Abort("invalid MPI placement, tag range or message size");
		// All ranks must agree on the mathematical plan and placement before kernels.
		std::vector<uint64_t> config={plan.Id(),plan.Shape().nx,plan.Shape().ny,plan.Shape().nz};
		for (const auto& d:plan.Domains()) {config.push_back(d.begin);config.push_back(d.end);}
		for (auto r:ranks) config.push_back(r);
		uint64_t n=config.size(),small=0,large=0;
		Check(MPI_Allreduce(&n,&small,1,MPI_UINT64_T,MPI_MIN,comm->value));
		Check(MPI_Allreduce(&n,&large,1,MPI_UINT64_T,MPI_MAX,comm->value));
		if (small!=large || n>INT_MAX) comm->Abort("ranks disagree on plan size");
		std::vector<uint64_t> lo(config.size()),hi(config.size());
		Check(MPI_Allreduce(config.data(),lo.data(),int(n),MPI_UINT64_T,MPI_MIN,comm->value));
		Check(MPI_Allreduce(config.data(),hi.data(),int(n),MPI_UINT64_T,MPI_MAX,comm->value));
		if (lo!=hi) comm->Abort("ranks disagree on plan or placement");
	}
	const HaloRecipe& Recipe(const MessageId& id,bool receive)
	{
		MainThread();
		if (id.epoch.plan!=plan.Id() || id.channel>=plan.Halos().size())
			throw std::logic_error("message outside MPI partition plan");
		const auto& r=plan.Halos()[id.channel];
		if (r.source!=id.source || r.destination!=id.destination || r.phase!=id.epoch.phase ||
			ranks[receive?id.destination:id.source]!=uint32_t(rank))
			throw std::logic_error("wrong MPI message endpoint or phase");
		return r;
	}
};

MPIHaloTransport::MPIHaloTransport(MPI_Comm source,const PartitionPlan& plan,const std::vector<uint32_t>& ranks)
	:m_Impl(std::make_shared<Impl>(source,plan,ranks)) {}
MPIHaloTransport::~MPIHaloTransport() {}

CompletionPtr MPIHaloTransport::Receive(const MessageId& id,HostBuffer destination)
{
	const auto& recipe=m_Impl->Recipe(id,true);
	if (!destination || destination->size()!=recipe.Samples()) throw std::invalid_argument("wrong MPI receive buffer");
	auto token=std::make_shared<MPICompletion>(m_Impl->comm);
	token->expected=MakeHeader(id,destination->size());token->header.fill(0);
	token->receive=true;token->buffer=destination;
	const int peer=m_Impl->ranks[id.source],tag=int(id.channel*2);
	Check(MPI_Irecv(token->header.data(),9,MPI_UINT64_T,peer,tag,m_Impl->comm->value,&token->requests[0]));
	Check(MPI_Irecv(destination->data(),int(destination->size()),MPI_FLOAT,peer,tag+1,m_Impl->comm->value,&token->requests[1]));
	m_Impl->pending.push_back(token);return token;
}
CompletionPtr MPIHaloTransport::Send(const MessageId& id,HostBuffer source)
{
	const auto& recipe=m_Impl->Recipe(id,false);
	if (!source || source->size()!=recipe.Samples()) throw std::invalid_argument("wrong MPI send buffer");
	auto token=std::make_shared<MPICompletion>(m_Impl->comm);
	token->header=MakeHeader(id,source->size());token->buffer=source;
	const int peer=m_Impl->ranks[id.destination],tag=int(id.channel*2);
	Check(MPI_Isend(token->header.data(),9,MPI_UINT64_T,peer,tag,m_Impl->comm->value,&token->requests[0]));
	Check(MPI_Isend(source->data(),int(source->size()),MPI_FLOAT,peer,tag+1,m_Impl->comm->value,&token->requests[1]));
	m_Impl->pending.push_back(token);return token;
}
void MPIHaloTransport::Progress()
{
	for (auto it=m_Impl->pending.begin();it!=m_Impl->pending.end();)
	{
		auto token=it->lock();
		if (!token || token->Ready()) it=m_Impl->pending.erase(it); else ++it;
	}
}
void MPIHaloTransport::Fail(const std::string& reason) {m_Impl->comm->Abort(reason);}

MPICollectiveContext::MPICollectiveContext(MPI_Comm source):m_Comm(MPI_COMM_NULL)
{MainThread();Check(MPI_Comm_dup(source,&m_Comm));Check(MPI_Comm_set_errhandler(m_Comm,MPI_ERRORS_RETURN));}
MPICollectiveContext::~MPICollectiveContext()
{int finalized=0;MPI_Finalized(&finalized);if(!finalized && m_Comm!=MPI_COMM_NULL) MPI_Comm_free(&m_Comm);}
double MPICollectiveContext::Reduce(double value,MPI_Op operation)
{
	if (!All(std::isfinite(value))) Fail("nonfinite global reduction input");
	double output=0;Check(MPI_Allreduce(&value,&output,1,MPI_DOUBLE,operation,m_Comm));
	if (!std::isfinite(output)) Fail("nonfinite global reduction result");
	return output;
}
double MPICollectiveContext::Minimum(double value) {return Reduce(value,MPI_MIN);}
double MPICollectiveContext::Sum(double value) {return Reduce(value,MPI_SUM);}
bool MPICollectiveContext::All(bool value)
{MainThread();int input=value?1:0,output=0;Check(MPI_Allreduce(&input,&output,1,MPI_INT,MPI_MIN,m_Comm));return output!=0;}
void MPICollectiveContext::Fail(const std::string& reason)
{std::cerr<<"Distributed run failed: "<<reason<<std::endl;MPI_Abort(m_Comm,72);std::terminate();}

} }
