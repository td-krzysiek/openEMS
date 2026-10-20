// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_DISTRIBUTED_TEST_FIXTURE_H
#define OPENEMS_DISTRIBUTED_TEST_FIXTURE_H
#include "FDTD/distributed/cpu_executor.h"
#include <cmath>
#include <stdexcept>

namespace test_fixture {
using namespace openems::distributed;

inline LocalOperatorData Data(const Domain& d)
{
	LocalOperatorData data;data.domain=d;data.timestep=1e-12;
	for (auto& a:data.coefficients) a.resize(3*d.OwnedNodes());
	for (auto& a:data.initial_fields) a.resize(3*d.OwnedNodes());
	const auto& g=d.grid;const size_t o=d.OwnedNodes();
	for (uint32_t c=0;c<3;++c)
		for (uint32_t x=d.begin;x<d.end;++x)
			for (uint32_t y=0;y<g.ny;++y)
				for (uint32_t z=0;z<g.nz;++z)
				{
					const size_t j=c*o+size_t(x-d.begin)*g.Plane()+size_t(y)*g.nz+z;
					const bool wall=x==0||x==g.nx-1||y==0||y==g.ny-1||z==0||z==g.nz-1;
					const bool pec=x==g.nx/3 && y>g.ny/2 && z>g.nz/2;
					data.coefficients[0][j]=(wall||pec)?0:1;
					data.coefficients[1][j]=(wall||pec)?0:((x/3+y/2+z/4)%3?.16f:.12f);
					data.coefficients[2][j]=1;
					data.coefficients[3][j]=.2f;
					const double dx=(double(x)-g.nx*.47)/3.1,dy=(double(y)-g.ny*.45)/2.1,dz=(double(z)-g.nz*.43)/1.7;
					const float pulse=float(std::exp(-dx*dx-dy*dy-dz*dz));
					data.initial_fields[0][j]=(wall||pec)?0:pulse*float(c+1)*.3f;
					data.initial_fields[1][j]=(wall||pec)?0:pulse*float(3-c)*.07f;
				}
	return data;
}

inline std::vector<Placement> Placements(const PartitionPlan& p,bool one_rank_per_domain=false)
{
	std::vector<Placement> result;
	for (const auto& d:p.Domains())
		result.push_back({d.id,one_rank_per_domain?d.id:0,"test-host","test","resource-"+std::to_string(d.id),size_t(8)<<30,0});
	return result;
}

inline double Error(const PartitionPlan& plan,const CPUDomainExecutor& reference,
	const std::vector<std::array<std::vector<float>,2>>& fields)
{
	double error=0,peak=0;
	for (size_t f=0;f<2;++f)
	{
		const auto full=reference.OwnedField(f?FieldPhase::Current:FieldPhase::Voltage);
		for (const auto& d:plan.Domains())
		{
			const auto& local=fields[d.id][f];
			if (local.size()!=3*d.OwnedNodes()) throw std::runtime_error("wrong owned-field output size");
			for (uint32_t c=0;c<3;++c)
				for (size_t j=0;j<d.OwnedNodes();++j)
				{
					const float a=full[c*d.grid.Nodes()+size_t(d.begin)*d.grid.Plane()+j];
					const float b=local[c*d.OwnedNodes()+j];
					if (!std::isfinite(a)||!std::isfinite(b)) throw std::runtime_error("nonfinite field");
					peak=std::max(peak,std::abs(double(a)));error=std::max(error,std::abs(double(a)-b));
				}
		}
	}
	if (peak<1e-15) throw std::runtime_error("trivial zero reference");
	return error/peak;
}
}
#endif
