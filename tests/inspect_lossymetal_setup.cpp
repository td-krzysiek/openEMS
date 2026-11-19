// SPDX-License-Identifier: GPL-3.0-or-later
#include <ContinuousStructure.h>
#include <tinyxml.h>
#include <FDTD/operator.h>
#include <FDTD/extensions/operator_ext_lossymetal.h>
#include <fstream>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <string>

// Reads only the protected state of its own extension instance.
class InspectLossyMetal : public Operator_Ext_LossyMetal
{
public:
    InspectLossyMetal(Operator* op, double fmax, bool& captured, unsigned int& surfaces)
        : Operator_Ext_LossyMetal(op, fmax), m_Captured(captured), m_Surfaces(surfaces)
    {
    }

    bool BuildExtension() override
    {
        const bool active = Operator_Ext_LossyMetal::BuildExtension();
        // CalcECOperator deletes inactive extensions. Capture while this
        // instance is still alive, including valid zero-surface fallbacks.
        Dump("ade.bin");
        m_Surfaces = m_numSurface;
        m_Captured = true;
        return active;
    }

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

private:
    bool& m_Captured;
    unsigned int& m_Surfaces;
};

int main(int argc, char** argv)
{
    if (argc < 3 || argc > 5) {
        std::cerr << "Usage: inspect_setup model.xml threads [--expect-surface] [--exclude-all]\n";
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
    bool expectSurface = false;
    bool excludeAll = false;
    for (int i = 3; i < argc; ++i) {
        const std::string option(argv[i]);
        if (option == "--expect-surface") {
            expectSurface = true;
        } else if (option == "--exclude-all") {
            excludeAll = true;
        } else {
            return 2;
        }
    }
    bool captured = false;
    unsigned int surfaces = 0;
    Excitation excitation;
    excitation.SetupGaussianPulse(20e9, 4e9);
    std::unique_ptr<Operator> op(Operator::New(std::stoi(argv[2])));
    op->SetEngineType(Engine_Types::BASIC);
    op->SetQuarterCellMaterialAvg();
    op->SetExcitationSignal(&excitation);
    op->SetTimeStepMethod(3);
    if (!op->SetGeometryCSX(&geometry)) {
        return 5;
    }
    int boundaries[6] = {0, 0, 0, 0, 0, 0};
    const char* names[6] = {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"};
    auto* fdtd = root->FirstChildElement("FDTD");
    if (fdtd == nullptr) {
        return 4;
    }
    auto* conditions = fdtd->FirstChildElement("BoundaryCond");
    for (int side = 0; side < 6; ++side) {
        const char* value = conditions ? conditions->Attribute(names[side]) : nullptr;
        if (value != nullptr && std::string(value) == "PMC") {
            boundaries[side] = 1;
        } else if (value != nullptr && std::string(value) != "PEC") {
            std::cerr << "Inspector supports PEC/PMC boundaries only\n";
            return 7;
        }
        // The regular openEMS orchestrator sets these explicitly too. The
        // low-level Operator API does not initialize boundary-layer sizes.
        op->SetBCSize(side, excludeAll ? op->GetNumberOfLines(side / 2, true) : 0);
    }
    op->SetBoundaryCondition(boundaries);
    auto* metal = new InspectLossyMetal(op.get(), excitation.GetMaxFreq(), captured, surfaces);
    op->AddExtension(metal);
    if (op->CalcECOperator(Operator::debugPEC | Operator::debugOperator) != 0) {
        return 6;
    }
    if (!captured || (expectSurface && surfaces == 0)) {
        std::cerr << "Missing expected LossyMetal surface state\n";
        return 8;
    }
    std::cout << "Captured LossyMetal surface edges: " << surfaces << '\n';
    std::cout << "ADE_INSPECTION_COMPLETE_NO_FDTD\n";
    return 0;
}
