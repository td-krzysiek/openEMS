// SPDX-License-Identifier: GPL-3.0-or-later
#include <ContinuousStructure.h>
#include <tinyxml.h>
#include <FDTD/operator.h>
#include <FDTD/extensions/operator_ext_lossymetal.h>
#include <fstream>
#include <memory>
#include <iostream>

// Reads only the protected state of its own extension instance.
class InspectLossyMetal : public Operator_Ext_LossyMetal
{
public:
    using Operator_Ext_LossyMetal::Operator_Ext_LossyMetal;
    void Dump(const char* filename) const
    {
        std::ofstream out(filename, std::ios::binary);
        auto write = [&](const auto& value) {
            out.write(reinterpret_cast<const char*>(&value), sizeof(value));
        };
        write(m_Order);
        write(m_numSurface);
        write(m_numFinite);
        write(m_numPEC);
        write(m_numH);
        write(m_numNoThickness);
        for (int o = 0; o < m_Order; ++o) {
            write(m_LM_Count[o]);
            write(m_volt_ADE_On[o]);
            write(m_curr_ADE_On[o]);
            write(m_volt_Lor_ADE_On[o]);
            write(m_curr_Lor_ADE_On[o]);
            for (int n = 0; n < 3; ++n) {
                for (size_t i = 0; i < m_LM_Count[o]; ++i) {
                    write(m_LM_pos[o][n][i]);
                    write(v_int_ADE[o][n][i]);
                    write(v_ext_ADE[o][n][i]);
                }
            }
        }
        for (const auto& fit : m_Fits) {
            write(fit.G);
            for (double value : fit.a) {
                write(value);
            }
        }
        if (!out) {
            throw std::runtime_error("Cannot write ADE diagnostic");
        }
    }
};

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: inspect_setup model.xml threads\n";
        return 2;
    }
    TiXmlDocument document(argv[1]);
    if (!document.LoadFile()) {
        return 3;
    }
    auto* root = document.FirstChildElement("openEMS");
    if (root == nullptr) {
        return 4;
    }
    ContinuousStructure geometry;
    const char* message = geometry.ReadFromXML(root);
    if ((message && *message) || !geometry.Update().empty()) {
        return 4;
    }
    Excitation excitation;
    excitation.SetupGaussianPulse(24.125e9, 4e9);
    std::unique_ptr<Operator> op(Operator::New(std::stoi(argv[2])));
    op->SetEngineType(Engine_Types::BASIC);
    op->SetQuarterCellMaterialAvg();
    op->SetExcitationSignal(&excitation);
    op->SetTimeStepMethod(3);
    if (!op->SetGeometryCSX(&geometry)) {
        return 5;
    }
    auto* metal = new InspectLossyMetal(op.get(), excitation.GetMaxFreq());
    op->AddExtension(metal);
    if (op->CalcECOperator(Operator::debugPEC | Operator::debugOperator) != 0) {
        return 6;
    }
    metal->Dump("ade.bin");
    std::cout << "ADE_INSPECTION_COMPLETE_NO_FDTD\n";
    return 0;
}
