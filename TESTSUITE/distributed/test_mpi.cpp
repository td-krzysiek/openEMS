// SPDX-License-Identifier: GPL-3.0-or-later
#include "fixture.h"
#include "FDTD/distributed/mpi_transport.h"
#include <iostream>
#include <memory>

using namespace openems::distributed;
int main(int argc,char** argv)
{
	int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
	int rank=0,size=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
	try
	{
		{
			const GridShape shape={37,21,14};
			// Two domains per rank exercise local and remote routes in one plan.
			std::vector<uint32_t> cuts,ranks;
			for(int n=0;n<=size*2;++n) cuts.push_back(uint32_t(uint64_t(shape.nx)*n/(size*2)));
			for(int n=0;n<size*2;++n) ranks.push_back(uint32_t(n/2));
			PartitionPlan p(7,shape,cuts);
			MPIHaloTransport transport(MPI_COMM_WORLD,p,ranks);MPICollectiveContext collectives(MPI_COMM_WORLD);
			std::vector<std::unique_ptr<CPUDomainExecutor>> owners;std::vector<DomainExecutor*> local;
			for (const auto& d:p.Domains()) if (ranks[d.id]==uint32_t(rank))
			{owners.emplace_back(new CPUDomainExecutor(test_fixture::Data(d)));local.push_back(owners.back().get());}
			FDTDCoordinator coordinator(p,transport,collectives,local);
			auto placement=test_fixture::Placements(p);for(auto& row:placement)row.process=ranks[row.domain];
			RunManager manager(p,placement,coordinator);manager.Run(120);
			PartitionPlan whole(1,shape,{0,shape.nx});CPUDomainExecutor reference(test_fixture::Data(whole.Domains()[0]));
			LocalHaloTransport plain;LocalCollectiveContext one;FDTDCoordinator ref(whole,plain,one,{&reference});
			RunManager reference_manager(whole,test_fixture::Placements(whole),ref);reference_manager.Run(120);
			double local_error=0,local_peak=0;
			for(const auto& executor:owners)
				for(int f=0;f<2;++f)
				{
					const auto phase=f?FieldPhase::Current:FieldPhase::Voltage;
					const auto a=executor->OwnedField(phase),b=reference.OwnedField(phase);
					const auto& d=p.Domains()[executor->Id()];
					for(uint32_t c=0;c<3;++c)for(size_t j=0;j<d.OwnedNodes();++j)
					{
						const float expected=b[c*shape.Nodes()+size_t(d.begin)*shape.Plane()+j];
						if(!std::isfinite(a[c*d.OwnedNodes()+j]))throw std::runtime_error("nonfinite MPI field");
						local_error=std::max(local_error,std::abs(double(a[c*d.OwnedNodes()+j])-expected));
						local_peak=std::max(local_peak,std::abs(double(expected)));
					}
				}
			if(!collectives.All(local_error==0 && local_peak>1e-15))throw std::runtime_error("MPI split mismatch");
			if(collectives.Minimum(rank+1.)!=1 || collectives.Sum(rank+1.)!=double(size*(size+1))/2)
				throw std::runtime_error("MPI reduction mismatch");
			if(argc>1 && std::string(argv[1])=="--fail-epoch")
			{
				const auto& recipe=p.Halos()[2];
				auto buffer=std::make_shared<std::vector<float>>(recipe.Samples());
				MessageId id={{7,uint64_t(rank==0?9:8),FieldPhase::Voltage},recipe.channel,recipe.source,recipe.destination};
				if(rank==int(ranks[recipe.destination]))transport.Receive(id,buffer)->Wait();
				if(rank==int(ranks[recipe.source]))transport.Send(id,buffer)->Wait();
				// An unexpectedly accepted epoch returns success, so WILL_FAIL fails.
				if(rank==0)std::cout<<"ERROR: wrong-epoch control was accepted\n";
			}
			if(rank==0)std::cout<<"PASS: "<<size<<" MPI ranks, "<<p.Domains().size()<<" domains, exact six-field agreement\n";
		}
		MPI_Finalize();return 0;
	}
	catch(const std::exception& e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,73);return 73;}
}
