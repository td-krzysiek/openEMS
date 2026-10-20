// SPDX-License-Identifier: GPL-3.0-or-later
#include "operator_fixture.h"
#include "FDTD/distributed/mpi_transport.h"
#ifdef OPENEMS_DOMAIN_TEST_HIP
#include "FDTD/hip/hip_domain_executor.h"
#endif
#include <iostream>
#include <memory>

using namespace openems::distributed;
int main(int argc,char** argv)
{
	int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
	int rank=0,size=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
	try
	{
		if(size!=2||argc<2||argc>3)throw std::invalid_argument("two ranks and a coefficient fixture are required");
		{
			const bool mixed=argc==3 && std::string(argv[2])=="--cpu-gpu";
			const auto full=test_fixture::LoadOperator(argv[1]);const auto shape=full.domain.grid;
			PartitionPlan p(9,shape,{0,shape.nx/3,shape.nx});
			MPIHaloTransport transport(MPI_COMM_WORLD,p,{0,1});MPICollectiveContext collectives(MPI_COMM_WORLD);
			const auto data=test_fixture::Slice(full,p.Domains()[rank]);
			std::unique_ptr<CPUDomainExecutor> cpu;
#ifdef OPENEMS_DOMAIN_TEST_HIP
			std::unique_ptr<HIPDomainExecutor> gpu;
#endif
			DomainExecutor* executor;
			if(mixed && rank==0){cpu.reset(new CPUDomainExecutor(data));executor=cpu.get();}
			else
			{
#ifdef OPENEMS_DOMAIN_TEST_HIP
				gpu.reset(new HIPDomainExecutor(data,rank));executor=gpu.get();
#else
				throw std::runtime_error("CPU-only executable assigned a GPU domain");
#endif
			}
			FDTDCoordinator coordinator(p,transport,collectives,{executor});
			RunManager manager(p,test_fixture::Placements(p,true),coordinator);manager.Run(120);
			PartitionPlan whole(1,shape,{0,shape.nx});CPUDomainExecutor reference(full);
			LocalHaloTransport rt;LocalCollectiveContext rc;FDTDCoordinator reference_coordinator(whole,rt,rc,{&reference});
			RunManager reference_manager(whole,test_fixture::Placements(whole),reference_coordinator);reference_manager.Run(120);
			double local_error=0,local_peak=0;
			const auto& d=p.Domains()[rank];
			for(int f=0;f<2;++f)
			{
				const auto phase=f?FieldPhase::Current:FieldPhase::Voltage;
				std::vector<float> values;
#ifdef OPENEMS_DOMAIN_TEST_HIP
				if(gpu)values=gpu->OwnedField(phase);else
#endif
				values=cpu->OwnedField(phase);
				const auto expected=reference.OwnedField(phase);
				for(uint32_t c=0;c<3;++c)for(size_t j=0;j<d.OwnedNodes();++j)
				{
					const float a=values[c*d.OwnedNodes()+j],b=expected[c*shape.Nodes()+size_t(d.begin)*shape.Plane()+j];
					if(!std::isfinite(a)||!std::isfinite(b))throw std::runtime_error("nonfinite MPI/HIP field");
					local_error=std::max(local_error,std::abs(double(a)-b));local_peak=std::max(local_peak,std::abs(double(b)));
				}
			}
			double error=0,peak=0;MPI_Allreduce(&local_error,&error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			MPI_Allreduce(&local_peak,&peak,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
			if(peak<1e-15||error/peak>2e-5)throw std::runtime_error("MPI/HIP field mismatch");
			if(rank==0)std::cout<<"PASS: two MPI ranks, "<<(mixed?"CPU+GPU":"GPU+GPU")
				<<", real operator, normalized six-field error="<<error/peak<<'\n';
		}
		MPI_Finalize();return 0;
	}
	catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,74);return 74;}
}
