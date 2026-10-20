// SPDX-License-Identifier: GPL-3.0-or-later
#include "FDTD/distributed/runtime.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace openems::distributed;
static void Require(bool yes, const char* message)
{if (!yes) throw std::runtime_error(message);}
template<typename F> static void Reject(F f)
{
	bool threw=false; try {f();} catch (const std::exception&) {threw=true;}
	Require(threw,"invalid input did not fail");
}

static void Plans()
{
	for (uint32_t nx=2;nx<38;++nx)
		for (uint32_t cut=1;cut<nx;++cut)
		{
			PartitionPlan p(12,{nx,7,5},{0,cut,nx});
			std::vector<unsigned int> owners(nx,0);
			for (const auto& d:p.Domains())
				for (uint32_t x=d.begin;x<d.end;++x)
				{
					Require(p.Owner(x)==d.id,"wrong sample owner"); ++owners[x];
					Require(d.ToLocal(x)==x-d.begin+1,"wrong local sample map");
					for (uint32_t c=0;c<3;++c)
						Require(d.LocalIndex(c,d.ToLocal(x),6,4)<3*d.StorageNodes(),"field index overflow");
				}
			Require(std::all_of(owners.begin(),owners.end(),[](unsigned int n){return n==1;}),"ownership gap/overlap");
			const auto& v=p.Halos()[0]; const auto& i=p.Halos()[1];
			Require(v.source==1 && v.destination==0 && v.source_x==1 && v.destination_x==cut+1,"wrong V halo direction");
			Require(i.source==0 && i.destination==1 && i.source_x==cut && i.destination_x==0,"wrong I halo direction");
			Require(v.Samples()==70 && i.Samples()==70,"wire component count changed");
		}
	PartitionPlan many(1,{13,5,3},{0,1,4,7,13});
	Require(many.Halos().size()==6,"multi-domain neighbour recipes missing");
	Reject([](){PartitionPlan p(0,{5,4,3},{0,5});});
	Reject([](){PartitionPlan p(1,{5,4,3},{0,2,2,5});});
	Reject([](){PartitionPlan p(1,{5,4,3},{0,3,2,5});});
	Reject([](){PartitionPlan p(1,{5,4,3},{1,5});});
	Reject([](){PartitionPlan p(1,{5,1,3},{0,5});});
	Reject([](){GridShape{UINT32_MAX,UINT32_MAX,UINT32_MAX}.Validate();});
	Reject([](){CheckedAdd(std::numeric_limits<size_t>::max(),1);});
	Reject([](){CheckedMultiply(std::numeric_limits<size_t>::max(),2);});
	Reject([&](){many.Owner(13);});
	Reject([&](){many.Domains()[1].ToLocal(0);});
}

static void Memory()
{
	PartitionPlan p(7,{10,5,3},{0,4,10});
	size_t b=BasicDomainBytes(p.Domains()[0])+BasicDomainBytes(p.Domains()[1]);
	std::vector<Placement> placement={{0,0,"host","cpu","shared-ram",b,0},
		{1,1,"host","cpu","shared-ram",b,0}};
	ValidatePlacement(p,placement);
	auto invalid=placement; invalid[1].domain=0;
	Reject([&](){ValidatePlacement(p,invalid);});
	invalid=placement; for(auto& row:invalid) row.memory_budget_bytes=b-1;
	Reject([&](){ValidatePlacement(p,invalid);});
	invalid=placement; invalid[0].auxiliary_bytes=1;
	Reject([&](){ValidatePlacement(p,invalid);});
	invalid=placement; invalid[1].backend="metal"; invalid[1].auxiliary_bytes=1;
	Reject([&](){ValidatePlacement(p,invalid);}); // shared physical RAM, different backend
}

static MessageId Id(uint64_t step=0)
{return {{2,step,FieldPhase::Current},1,0,1};}

static void Transfers()
{
	CompletionPtr send,receive;
	auto source=std::make_shared<std::vector<float>>(4,3.5f);
	auto destination=std::make_shared<std::vector<float>>(4,-9.f);
	std::weak_ptr<std::vector<float>> lifetime=source;
	{
		LocalHaloTransport t(5);
		receive=t.Receive(Id(),destination);
		send=t.Send(Id(),source); source.reset();
		Require(!lifetime.expired(),"send buffer released before completion");
		Require(!receive->Ready() && destination->at(0)==-9.f,"delayed receive completed too early");
		receive->Wait(); send->Wait();
		Require(destination->at(3)==3.5f && lifetime.expired(),"copy or send release lifetime failed");
	}
	Require(send->Ready() && receive->Ready(),"completed tokens require live transport");
	{
		LocalHaloTransport t;
		auto a=t.Receive(Id(),destination);
		Reject([&](){t.Receive(Id(),destination);});
		t.Fail("injected copy failure"); Reject([&](){a->Wait();});
	}
	{
		LocalHaloTransport t;
		auto a=t.Receive(Id(),destination);
		auto b=t.Send(Id(),std::make_shared<std::vector<float>>(3,0));
		Reject([&](){a->Wait();});
		t.Fail("mismatch stopped"); Reject([&](){b->Wait();});
	}
	{
		LocalHaloTransport t(20);
		auto a=t.Receive(Id(0),destination);
		auto b=t.Send(Id(1),destination);
		Require(!a->Ready() && !b->Ready(),"wrong epoch matched");
		t.Fail("wrong epoch deliberately cancelled");
		Reject([&](){a->Wait();}); Reject([&](){b->Wait();});
	}
}

class DelayedUnpack : public Completion
{
public:
	explicit DelayedUnpack(std::function<void()> action):m_Action(action),m_Polls(0),m_Done(false) {}
	bool Ready() override {if(!m_Done && ++m_Polls>3){m_Action();m_Done=true;}return m_Done;}
	void Wait() override {while(!Ready()) {}}
private:
	std::function<void()> m_Action;int m_Polls;bool m_Done;
};

struct MockExecutor : DomainExecutor
{
	uint32_t id; std::vector<std::string>& trace;
	bool fail=false; double dt=1.;
	MockExecutor(uint32_t n,std::vector<std::string>& t):id(n),trace(t) {}
	uint32_t Id() const override {return id;}
	double Timestep() const override {return dt;}
	CompletionPtr Update(FieldPhase phase,uint64_t step) override
	{
		if (fail) throw std::runtime_error("injected executor failure");
		trace.push_back(std::string(phase==FieldPhase::Voltage?"V":"I")+std::to_string(id)+":"+std::to_string(step));
		return Completed();
	}
	CompletionPtr Pack(const HaloRecipe& r,HostBuffer b) override
	{Require(r.source==id,"pack not owned"); std::fill(b->begin(),b->end(),float(id));return Completed();}
	CompletionPtr Unpack(const HaloRecipe& r,HostBuffer b) override
	{return std::make_shared<DelayedUnpack>([this,r,b](){Require(r.destination==id && b->at(0)==float(r.source),"wrong received owner/sample");trace.push_back("unpack");});}
	void Observe(FieldPhase phase,uint64_t) override
	{trace.push_back(phase==FieldPhase::Voltage?"observeV":"observeI");}
};

static void Coordination()
{
	PartitionPlan p(1,{9,5,3},{0,3,9});
	LocalHaloTransport transport(10); LocalCollectiveContext collectives;
	std::vector<std::string> trace; MockExecutor a(0,trace),b(1,trace);
	FDTDCoordinator coordinator(p,transport,collectives,{&a,&b});
	Reject([&](){coordinator.Step();});
	std::vector<Placement> placement={{0,0,"local","test","0",1000000,0},{1,0,"local","test","1",1000000,0}};
	RunManager manager(p,placement,coordinator);manager.Run(2);
	Require(manager.Finished() && coordinator.NextStep()==2,"global timestep did not advance exactly once");
	std::vector<std::string> expected={"unpack",
		"V0:0","V1:0","unpack","observeV","observeV","I0:0","I1:0","unpack","observeI","observeI",
		"V0:1","V1:1","unpack","observeV","observeV","I0:1","I1:1","unpack","observeI","observeI"};
	Require(trace==expected,"half-step / receive / observation ordering changed");
	Reject([&](){manager.Run(1);});
	Reject([&](){coordinator.Bootstrap();});
	HaloExchange x(p,transport,{&a,&b});
	Reject([&](){x.Begin({999,0,FieldPhase::Current});});
	Reject([&](){x.Begin({1,0,FieldPhase::Voltage});});
	auto initial=x.Begin({1,0,FieldPhase::Current});
	Reject([&](){x.Begin({1,1,FieldPhase::Voltage});}); x.Wait(initial);
	Reject([&](){x.Begin({1,0,FieldPhase::Current});});
	Reject([&](){x.Begin({1,2,FieldPhase::Voltage});});
	LocalHaloTransport failed_transport;
	FDTDCoordinator failed(p,failed_transport,collectives,{&a,&b});
	RunManager failed_manager(p,placement,failed);failed_manager.Start();a.fail=true;
	Reject([&](){failed_manager.Advance(1);});Require(failed.NextStep()==0,"failed epoch advanced global step");
	Reject([&](){failed.Step();});
	Reject([&](){failed_manager.Finish();});
	a.fail=false;
	LocalHaloTransport missing_transport;
	FDTDCoordinator missing(p,missing_transport,collectives,{&a});
	Reject([&](){missing.Bootstrap();});
	Require(missing.Failed(),"missing executor did not fail the coordinator");
	Reject([&](){FDTDCoordinator duplicate(p,missing_transport,collectives,{&a,&a,&b});});
	LocalHaloTransport timestep_transport;b.dt=2.;
	FDTDCoordinator timestep(p,timestep_transport,collectives,{&a,&b});
	Reject([&](){timestep.Bootstrap();});
	Require(timestep.Failed(),"different global timesteps were accepted");
	b.dt=1.;
	Reject([&](){collectives.Sum(std::numeric_limits<double>::quiet_NaN());});
}

int main()
{
	try {Plans();Memory();Transfers();Coordination();}
	catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
	std::cout<<"PASS: ownership, memory, delayed transport, buffer lifetime, failure and half-step scheduling\n";
}
