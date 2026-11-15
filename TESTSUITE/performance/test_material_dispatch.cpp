// SPDX-License-Identifier: GPL-3.0-or-later
// Downstream subclasses overriding only the legacy hook must still be called.
#ifdef NDEBUG
#error "Compile with assertions enabled"
#endif
#include "FDTD/extensions/engine_ext_lorentzmaterial.h"
#include "FDTD/extensions/operator_ext_lorentzmaterial.h"
#include <cassert>
#include <iostream>

class CustomMaterial : public Engine_Ext_LorentzMaterial
{
public:
    explicit CustomMaterial(Operator_Ext_LorentzMaterial* op)
        : Engine_Ext_LorentzMaterial(op), voltage_calls(0), current_calls(0) {}
    void DoPreVoltageUpdates() override {++voltage_calls;}
    void DoPreCurrentUpdates() override {++current_calls;}
    unsigned int voltage_calls, current_calls;
};

int main()
{
    Operator_Ext_LorentzMaterial op(NULL); // zero-order operator; no solver needed
    CustomMaterial custom(&op);
    custom.SetNumberOfThreads(16);
    Engine_Extension& extension = custom;
    for (int worker=0; worker<16; ++worker)
    {
        extension.DoPreVoltageUpdates(worker);
        extension.DoPreCurrentUpdates(worker);
    }
    assert(custom.voltage_calls==1 && custom.current_calls==1);
    std::cout << "Legacy material subclass dispatch checks passed\n";
}
