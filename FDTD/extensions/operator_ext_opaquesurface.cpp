// SPDX-License-Identifier: GPL-3.0-or-later
#include "operator_ext_opaquesurface.h"
#include "engine_ext_opaquesurface.h"
#include "FDTD/operator.h"
#include "CSPrimBox.h"
#include "CSPropMaterial.h"
#include "CSTransform.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace {
struct Box {double lo[3],hi[3]; unsigned int fit;};
bool interior(const double* c,const Box& b) {
    for (int n=0;n<3;++n) if (c[n]<=b.lo[n] || c[n]>=b.hi[n]) return false;
    return true;
}
double number(const std::string& s) {
    size_t used=0; double v=std::stod(s,&used);
    if (used!=s.size() || !std::isfinite(v)) throw std::runtime_error("Invalid opaque surface coefficient");
    return v;
}
}

void Operator_Ext_OpaqueSurface::PreprocessEC()
{
    fits.clear(); edges.clear(); interior_h.clear(); pml_guard_edges=0;
    auto* csx=m_Op->GetGeometryCSX();
    std::vector<Box> boxes;
    for (size_t i=0;i<csx->GetQtyProperties();++i) {
        auto* prop=csx->GetProperty(i);
        if (!prop->ExistAttribute("ExperimentalOpaqueSurfaceV1")) continue;
        if (!(prop->GetType() & CSProperties::METAL)) throw std::runtime_error("Opaque surface requires Metal boxes");
        Fit f;
        std::stringstream text(prop->GetAttributeValue("ExperimentalOpaqueSurfaceV1")); std::string token;
        if (!std::getline(text,token,';')) throw std::runtime_error("Missing opaque surface fit");
        f.g=number(token);
        while (std::getline(text,token,';')) {
            size_t comma=token.find(',');
            if (comma==std::string::npos) throw std::runtime_error("Invalid opaque surface branch");
            f.w.push_back(number(token.substr(0,comma))); f.tau.push_back(number(token.substr(comma+1)));
            if (f.w.back()<=0 || f.tau.back()<=0) throw std::runtime_error("Opaque surface branches must be positive");
        }
        if (f.g<0 || f.w.empty() || f.w.size()>128) throw std::runtime_error("Invalid opaque surface fit order/conductance");
        unsigned int fi=fits.size(); fits.push_back(f);
        for (auto* primitive:prop->GetAllPrimitives()) {
            auto* box=dynamic_cast<CSPrimBox*>(primitive);
            if (!box || primitive->GetTransform()->HasTransform()) throw std::runtime_error("Opaque surface supports untransformed Cartesian boxes only");
            Box b; b.fit=fi;
            for (int n=0;n<3;++n) {
                b.lo[n]=std::min(box->GetCoord(2*n),box->GetCoord(2*n+1));
                b.hi[n]=std::max(box->GetCoord(2*n),box->GetCoord(2*n+1));
                if (!(b.hi[n]>b.lo[n])) throw std::runtime_error("Opaque surface requires positive physical thickness");
                for (auto value : {b.lo[n],b.hi[n]}) {
                    bool inside=false; auto pos=m_Op->SnapToMeshLine(n,value,inside);
                    if (!inside || std::abs(m_Op->GetDiscLine(n,pos)-value)>1e-8)
                        throw std::runtime_error("Opaque surface box faces must coincide with mesh lines");
                }
            }
            for (const auto& other:boxes) {
                bool overlap=true;
                for (int n=0;n<3;++n) overlap=overlap && b.hi[n]>=other.lo[n] && other.hi[n]>=b.lo[n];
                if (overlap) throw std::runtime_error("Opaque surface boxes must be disjoint, including their faces");
            }
            boxes.push_back(b);
        }
    }
    if (boxes.empty()) return;
    if (m_Op_Cyl) throw std::runtime_error("Opaque surface requires Cartesian coordinates");
    const double unit=m_Op->GetGridDelta(), eps0=8.854187817620389e-12;
    unsigned int pos[3];
    // Exclude interior magnetic unknowns: otherwise finite face E drives a
    // spurious inductive field in what is supposed to be excluded metal.
    for (int n=0;n<3;++n)
    for (pos[0]=0;pos[0]+1<m_Op->GetNumberOfLines(0);++pos[0])
    for (pos[1]=0;pos[1]+1<m_Op->GetNumberOfLines(1);++pos[1])
    for (pos[2]=0;pos[2]+1<m_Op->GetNumberOfLines(2);++pos[2]) {
        double c[3]; m_Op->GetYeeCoords(n,pos,c,true);
        for (const auto& b:boxes) if (interior(c,b)) {
            interior_h.push_back({{static_cast<unsigned int>(n),pos[0],pos[1],pos[2]}});break;
        }
    }
    for (unsigned int n=0;n<3;++n)
    for (pos[0]=0;pos[0]<m_Op->GetNumberOfLines(0);++pos[0])
    for (pos[1]=0;pos[1]<m_Op->GetNumberOfLines(1);++pos[1])
    for (pos[2]=0;pos[2]<m_Op->GetNumberOfLines(2);++pos[2]) {
        if (pos[n]+1>=m_Op->GetNumberOfLines(n)) continue;
        double c[3]; m_Op->GetYeeCoords(n,pos,c,false);
        const unsigned int p=(n+1)%3,q=(n+2)%3;
        const double length=m_Op->GetEdgeLength(n,pos);
        double scale=0; unsigned int fi=0;
        for (const auto& b:boxes) {
            if (c[n]<=b.lo[n] || c[n]>=b.hi[n]) continue;
            for (unsigned int normal : {p,q}) {
                // A metal face on the outer domain wall has no exterior
                // half cell and is governed by the domain boundary condition.
                if (!pos[normal] || pos[normal]+1>=m_Op->GetNumberOfLines(normal)) continue;
                unsigned int along=3-n-normal;
                if (c[along]<b.lo[along]-1e-9 || c[along]>b.hi[along]+1e-9) continue;
                if (std::abs(c[normal]-b.lo[normal])>1e-9 && std::abs(c[normal]-b.hi[normal])>1e-9) continue;
                double low=pos[along] ? (m_Op->GetDiscLine(along,pos[along]-1)+c[along])/2 : c[along];
                double high=pos[along]+1<m_Op->GetNumberOfLines(along) ? (m_Op->GetDiscLine(along,pos[along]+1)+c[along])/2 : c[along];
                scale+=std::max(0.,std::min(high,b.hi[along])-std::max(low,b.lo[along]))*unit/length;
                fi=b.fit;
            }
        }
        if (scale==0) continue;
        bool pml=false;
        for (int axis=0;axis<3;++axis) {
            int low=m_Op->GetBCSize(2*axis),high=m_Op->GetBCSize(2*axis+1);
            pml=pml || (low>0 && pos[axis]<=static_cast<unsigned int>(low)) ||
                (high>0 && pos[axis]+high+1>=m_Op->GetNumberOfLines(axis));
        }
        if (pml) {++pml_guard_edges;continue;}
        // Integrate only exterior dielectric quadrants of the electric dual face.
        double cap=0,gd=0;
        for (int sp : {-1,1}) for (int sq : {-1,1}) {
            if ((sp<0 && !pos[p]) || (sq<0 && !pos[q]) ||
                (sp>0 && pos[p]+1>=m_Op->GetNumberOfLines(p)) || (sq>0 && pos[q]+1>=m_Op->GetNumberOfLines(q))) continue;
            double dp=std::abs(m_Op->GetDiscLine(p,pos[p]+sp)-c[p]);
            double dq=std::abs(m_Op->GetDiscLine(q,pos[q]+sq)-c[q]);
            double sample[3]={c[0],c[1],c[2]};sample[p]+=sp*dp/4;sample[q]+=sq*dq/4;
            bool metal=false;for (const auto& b:boxes) metal=metal || interior(sample,b);
            if (metal) continue;
            auto* mat=dynamic_cast<CSPropMaterial*>(csx->GetPropertyByCoordPriority(sample,CSProperties::MATERIAL));
            double eps=mat ? mat->GetEpsilonWeighted(n,sample):m_Op->GetBackgroundEpsR();
            double g=mat ? mat->GetKappaWeighted(n,sample):m_Op->GetBackgroundKappa();
            if (!std::isfinite(eps) || eps<=0 || !std::isfinite(g) || g<0) throw std::runtime_error("Opaque surface requires passive nondispersive exterior dielectric");
            cap+=eps0*eps*dp*dq*unit*unit/(4*length);
            gd+=g*dp*dq*unit*unit/(4*length);
        }
        if (!(cap>0)) throw std::runtime_error("Opaque surface has no exterior capacitance");
        Edge e;e.n=n;e.fit=fi;e.scale=scale;e.capacitance=cap;e.dielectric_g=gd;
        std::copy(pos,pos+3,e.pos);edges.push_back(e);
        // Use the same physical exterior C in the CFL calculation and update.
        // No artificial stabilizing capacitance is added to the RF model.
        const unsigned int index=m_Op->MainOp->SetPos(pos[0],pos[1],pos[2]);
        m_Op->EC_C[n][index]=cap;
        m_Op->EC_G[n][index]=gd;
    }
    if (edges.empty()) throw std::runtime_error("No active opaque surface edges");
}

bool Operator_Ext_OpaqueSurface::BuildExtension()
{
    if (fits.empty()) return false;
    dt=m_Op->GetTimestep();
    for (const auto& h:interior_h) {
        m_Op->SetII(h[0],h[1],h[2],h[3],0);
        m_Op->SetIV(h[0],h[1],h[2],h[3],0);
    }
    for (auto& e:edges) {
        const auto& fit=fits[e.fit];double bs=0;
        for (size_t j=0;j<fit.w.size();++j) bs+=dt*fit.w[j]*e.scale/(2*fit.tau[j]+dt);
        const double g=e.dielectric_g+fit.g*e.scale;
        const double denom=e.capacitance/dt+(g+bs)/2;
        e.vi=1/denom;e.vv=(e.capacitance/dt-(g+bs)/2)/denom;
        m_Op->SetVV(e.n,e.pos[0],e.pos[1],e.pos[2],0);
        m_Op->SetVI(e.n,e.pos[0],e.pos[1],e.pos[2],0);
    }
    std::cout << "EXPERIMENTAL_OPAQUE_SURFACE_V1 active: " << edges.size() << " edges, " << interior_h.size()
        << " excluded interior H, " << pml_guard_edges << " PML-guard edges kept PEC" << std::endl;
    return true;
}

Engine_Extension* Operator_Ext_OpaqueSurface::CreateEngineExtention()
{
    return new Engine_Ext_OpaqueSurface(this);
}
