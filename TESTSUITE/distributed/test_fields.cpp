// SPDX-License-Identifier: GPL-3.0-or-later
#include "fixture.h"
#include <iostream>
#include <memory>

using namespace openems::distributed;

class ZeroHalo : public LocalHaloTransport
{
public:
	CompletionPtr Send(const MessageId& id,HostBuffer source) override
	{return LocalHaloTransport::Send(id,std::make_shared<std::vector<float>>(source->size(),0));}
};

static double Solve(GridShape shape,std::vector<uint32_t> cuts,bool negative)
{
	PartitionPlan single(1,shape,{0,shape.nx}),split(2,shape,cuts);
	CPUDomainExecutor reference(test_fixture::Data(single.Domains()[0]));
	std::vector<std::unique_ptr<CPUDomainExecutor>> owned;
	std::vector<DomainExecutor*> executors;
	for (const auto& d:split.Domains())
	{owned.emplace_back(new CPUDomainExecutor(test_fixture::Data(d)));executors.push_back(owned.back().get());}
	LocalHaloTransport reference_transport,transport(17);ZeroHalo broken;
	LocalCollectiveContext collectives;
	FDTDCoordinator full(single,reference_transport,collectives,{&reference});
	FDTDCoordinator divided(split,negative?static_cast<HaloTransport&>(broken):transport,collectives,executors);
	RunManager a(single,test_fixture::Placements(single),full),b(split,test_fixture::Placements(split),divided);
	a.Run(120);b.Run(120);
	std::vector<std::array<std::vector<float>,2>> fields;
	for (const auto& e:owned) fields.push_back({{e->OwnedField(FieldPhase::Voltage),e->OwnedField(FieldPhase::Current)}});
	return test_fixture::Error(split,reference,fields);
}

int main()
{
	try
	{
		unsigned int cases=0;
		for (uint32_t nx:{35,36})
			for (uint32_t cut:{1u,11u,18u,25u,nx-1})
			{
				const double error=Solve({nx,19,13},{0,cut,nx},false);
				if (error!=0) throw std::runtime_error("CPU split fields differ");
				++cases;
			}
		if (Solve({37,21,14},{0,1,8,23,36,37},false)!=0) throw std::runtime_error("five-domain fields differ");
		const double negative=Solve({35,19,13},{0,17,35},true);
		if (negative<1e-3) throw std::runtime_error("missing-halo control was not detected");
		auto unsupported=test_fixture::Data(PartitionPlan(1,{5,5,5},{0,5}).Domains()[0]);
		unsupported.required_extensions.push_back("opaque");
		bool rejected=false;try {CPUDomainExecutor e(unsupported);}catch(const std::invalid_argument&){rejected=true;}
		if (!rejected) throw std::runtime_error("unsupported physics silently accepted");
		std::cout<<"PASS: "<<cases+1<<" CPU field cases; negative-control error="<<negative<<'\n';
	}
	catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
