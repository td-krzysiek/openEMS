// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPERATOR_EXT_OPAQUESURFACE_H
#define OPERATOR_EXT_OPAQUESURFACE_H
#include "operator_extension.h"
#include <array>
#include <vector>

// Experimental, grid-conforming Cartesian boxes only. Opaque RF limit.
class Operator_Ext_OpaqueSurface : public Operator_Extension
{
public:
    struct Fit { double g; std::vector<double> w, tau; };
    struct Edge {
        unsigned int n, pos[3], fit;
        double vv, vi, scale;
    };
    explicit Operator_Ext_OpaqueSurface(Operator* op) : Operator_Extension(op) {}
    virtual bool BuildExtension();
    virtual Engine_Extension* CreateEngineExtention();
    virtual std::string GetExtensionName() const {return "Experimental opaque surface V1";}
    std::vector<Fit> fits;
    std::vector<Edge> edges;
    double dt;
};
#endif
