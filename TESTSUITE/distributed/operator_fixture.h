// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_DOMAIN_OPERATOR_FIXTURE_H
#define OPENEMS_DOMAIN_OPERATOR_FIXTURE_H
#include "fixture.h"
#include <algorithm>
#include <fstream>

namespace test_fixture {
inline LocalOperatorData LoadOperator(const std::string& path)
{
	LocalOperatorData data;std::ifstream in(path);std::string magic;GridShape shape;
	in>>magic>>shape.nx>>shape.ny>>shape.nz>>data.timestep;
	if(!in || magic!="OPENEMS_DOMAIN_OPERATOR_V1")throw std::runtime_error("invalid coefficient fixture");
	shape.Validate();if(shape.Nodes()>2000000)throw std::runtime_error("coefficient fixture exceeds short-test bound");
	data.domain={0,0,shape.nx,shape};
	for(auto& a:data.coefficients){a.resize(3*shape.Nodes());for(auto& v:a)in>>v;}
	for(auto& a:data.initial_fields){a.resize(3*shape.Nodes());for(auto& v:a)in>>v;}
	if(!in)throw std::runtime_error("truncated coefficient fixture");
	in>>std::ws;if(!in.eof())throw std::runtime_error("trailing coefficient fixture data");
	data.Validate();return data;
}
inline LocalOperatorData Slice(const LocalOperatorData& full,const Domain& d)
{
	const auto& shape=full.domain.grid;
	if(shape.nx!=d.grid.nx||shape.ny!=d.grid.ny||shape.nz!=d.grid.nz)
		throw std::invalid_argument("fixture/partition grid mismatch");
	LocalOperatorData data;data.domain=d;data.timestep=full.timestep;
	const size_t plane=d.grid.Plane(),global=d.grid.Nodes(),owned=d.OwnedNodes();
	for(size_t n=0;n<6;++n)
	{
		auto& destination=n<4?data.coefficients[n]:data.initial_fields[n-4];
		const auto& source=n<4?full.coefficients[n]:full.initial_fields[n-4];
		destination.resize(3*owned);
		for(uint32_t c=0;c<3;++c)
			std::copy_n(source.data()+c*global+size_t(d.begin)*plane,owned,destination.data()+c*owned);
	}
	data.Validate();return data;
}
}
#endif
