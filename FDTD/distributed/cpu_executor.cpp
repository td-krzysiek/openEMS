// SPDX-License-Identifier: GPL-3.0-or-later
#include "cpu_executor.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace openems { namespace distributed {

void LocalOperatorData::Validate() const
{
	domain.grid.Validate();
	if (domain.begin>=domain.end || domain.end>domain.grid.nx || !std::isfinite(timestep) || timestep<=0)
		throw std::invalid_argument("invalid local operator or global timestep");
	if (!required_extensions.empty())
		throw std::invalid_argument("distributed Cartesian core does not yet support extension: "+required_extensions.front());
	const size_t count=CheckedMultiply(3,domain.OwnedNodes());
	for (const auto& a:coefficients)
		if (a.size()!=count || !std::all_of(a.begin(),a.end(),[](float v){return std::isfinite(v);}))
			throw std::invalid_argument("invalid final operator coefficient slice");
	for (const auto& a:initial_fields)
		if ((!a.empty() && a.size()!=count) || !std::all_of(a.begin(),a.end(),[](float v){return std::isfinite(v);}))
			throw std::invalid_argument("invalid initial owned fields");
}

CPUDomainExecutor::CPUDomainExecutor(const LocalOperatorData& data):m_Data(data)
{
	m_Data.Validate(); const auto& d=m_Data.domain;
	m_Voltage.resize(CheckedMultiply(3,d.StorageNodes()),0);
	m_Current.resize(m_Voltage.size(),0);
	for (uint32_t c=0;c<3;++c)
		for (uint32_t x=0;x<d.end-d.begin;++x)
			for (int field=0;field<2;++field)
				if (!data.initial_fields[field].empty())
					std::copy_n(data.initial_fields[field].data()+c*d.OwnedNodes()+x*d.grid.Plane(),d.grid.Plane(),
						(field?m_Current:m_Voltage).data()+d.LocalIndex(c,x+1,0,0));
	// Initial fields are not kept as another full replica after construction.
	for (auto& a:m_Data.initial_fields) std::vector<float>().swap(a);
}

void CPUDomainExecutor::Voltage()
{
	const auto& d=m_Data.domain; const auto& g=d.grid;
	const size_t s=d.StorageNodes(),p=g.Plane(),o=d.OwnedNodes();
	const auto& vv=m_Data.coefficients[0];const auto& vi=m_Data.coefficients[1];
	for (uint32_t x=1;x<=d.end-d.begin;++x)
		for (uint32_t y=0;y<g.ny;++y)
			for (uint32_t z=0;z<g.nz;++z)
			{
				const size_t i=size_t(x)*p+size_t(y)*g.nz+z,j=i-p;
				const size_t sx=(d.begin+x-1)?p:0,sy=y?g.nz:0,sz=z?1:0;
				m_Voltage[i]*=vv[j];
				m_Voltage[i]+=vi[j]*(m_Current[2*s+i]-m_Current[2*s+i-sy]-m_Current[s+i]+m_Current[s+i-sz]);
				m_Voltage[s+i]*=vv[o+j];
				m_Voltage[s+i]+=vi[o+j]*(m_Current[i]-m_Current[i-sz]-m_Current[2*s+i]+m_Current[2*s+i-sx]);
				m_Voltage[2*s+i]*=vv[2*o+j];
				m_Voltage[2*s+i]+=vi[2*o+j]*(m_Current[s+i]-m_Current[s+i-sx]-m_Current[i]+m_Current[i-sy]);
			}
}

void CPUDomainExecutor::Current()
{
	const auto& d=m_Data.domain;const auto& g=d.grid;
	const size_t s=d.StorageNodes(),p=g.Plane(),o=d.OwnedNodes();
	const auto& ii=m_Data.coefficients[2];const auto& iv=m_Data.coefficients[3];
	for (uint32_t x=1;x<=d.end-d.begin;++x)
	{
		if (d.begin+x-1==g.nx-1) continue;
		for (uint32_t y=0;y<g.ny-1;++y)
			for (uint32_t z=0;z<g.nz-1;++z)
			{
				const size_t i=size_t(x)*p+size_t(y)*g.nz+z,j=i-p;
				m_Current[i]*=ii[j];
				m_Current[i]+=iv[j]*(m_Voltage[2*s+i]-m_Voltage[2*s+i+g.nz]-m_Voltage[s+i]+m_Voltage[s+i+1]);
				m_Current[s+i]*=ii[o+j];
				m_Current[s+i]+=iv[o+j]*(m_Voltage[i]-m_Voltage[i+1]-m_Voltage[2*s+i]+m_Voltage[2*s+i+p]);
				m_Current[2*s+i]*=ii[2*o+j];
				m_Current[2*s+i]+=iv[2*o+j]*(m_Voltage[s+i]-m_Voltage[s+i+p]-m_Voltage[i]+m_Voltage[i+g.nz]);
			}
	}
}
CompletionPtr CPUDomainExecutor::Update(FieldPhase phase,uint64_t)
{if (phase==FieldPhase::Voltage) Voltage();else if(phase==FieldPhase::Current) Current();else throw std::invalid_argument("invalid field phase");return Completed();}
CompletionPtr CPUDomainExecutor::Pack(const HaloRecipe& r,HostBuffer b)
{
	const auto& d=m_Data.domain;
	if (!b || b->size()!=r.Samples() || r.source!=Id() || r.plane_samples!=d.grid.Plane() ||
		r.source_x==0 || r.source_x>d.end-d.begin) throw std::invalid_argument("invalid owned halo packing recipe");
	const auto& field=r.phase==FieldPhase::Voltage?m_Voltage:m_Current;
	for (uint32_t c=1;c<3;++c)
		std::copy_n(field.data()+d.LocalIndex(c,r.source_x,0,0),r.plane_samples,b->data()+(c-1)*r.plane_samples);
	return Completed();
}
CompletionPtr CPUDomainExecutor::Unpack(const HaloRecipe& r,HostBuffer b)
{
	const auto& d=m_Data.domain;
	if (!b || b->size()!=r.Samples() || r.destination!=Id() || r.plane_samples!=d.grid.Plane() ||
		(r.destination_x!=0 && r.destination_x!=d.end-d.begin+1)) throw std::invalid_argument("invalid ghost unpacking recipe");
	auto& field=r.phase==FieldPhase::Voltage?m_Voltage:m_Current;
	for (uint32_t c=1;c<3;++c)
		std::copy_n(b->data()+(c-1)*r.plane_samples,r.plane_samples,field.data()+d.LocalIndex(c,r.destination_x,0,0));
	return Completed();
}
std::vector<float> CPUDomainExecutor::OwnedField(FieldPhase phase) const
{
	const auto& d=m_Data.domain;const auto& a=phase==FieldPhase::Voltage?m_Voltage:m_Current;
	std::vector<float> out(3*d.OwnedNodes());
	for (uint32_t c=0;c<3;++c)
		std::copy_n(a.data()+d.LocalIndex(c,1,0,0),d.OwnedNodes(),out.data()+c*d.OwnedNodes());
	return out;
}

} }
