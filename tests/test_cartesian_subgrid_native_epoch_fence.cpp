/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "subgrid_native_gather_fixture.h"

static void Control(const SubgridDeviceTarget &target, unsigned int axis, bool loss, bool domains)
{
    const double unit = axis == 1 ? .0007 : .001, offset = axis == 2 ? -3 : 0;
    Fixture fenced(target, axis, loss, unit, offset, domains), control(target, axis, loss, unit, offset, domains);
    std::unique_ptr<Engine_Interface_FDTD> interfaces[2][2];
    std::unique_ptr<Engine_Field_Gather> gathers[2][2];
    int ids[2][2] = {};
    unsigned int count[3] = {3, 3, 3};
    unsigned int positions[3] = {6, 7, 8};
    unsigned int *lines[3] = {positions, positions, positions};
    Fixture *fixtures[2] = {&fenced, &control};
    for (unsigned int run = 0; run < 2; ++run)
    {
        for (unsigned int field = 0; field < 2; ++field)
        {
            interfaces[run][field] = std::make_unique<Engine_Interface_FDTD>(fixtures[run]->op.get());
            interfaces[run][field]->SetInterpolationType(Engine_Interface_Base::NODE_INTERPOLATE);
            gathers[run][field].reset(interfaces[run][field]->CreateFieldGather(field != 0, count, lines));
            ids[run][field] = interfaces[run][field]->CreateFieldDFT(gathers[run][field].get(), 2);
            Require(ids[run][field] >= 0, "native fence DFT preparation failed");
        }
        fixtures[run]->op->SealObservations();
    }
    Reject([&] { fenced.op->FenceObservationEpoch(1); });
    const std::size_t workers = domains ? 2 : 1;
    Require(fenced.op->FenceObservationEpoch(0) == workers, "startup native fence worker count differs");
    ArrayLib::ArrayNIJK<std::complex<float>> values[2][2];
    for (auto &run : values)
    {
        for (auto &field : run)
        {
            field.Init("FenceDFT", std::array<unsigned int, 3>{3, 3, 3});
        }
    }
    for (unsigned int batch = 0; batch < 8; ++batch)
    {
        fenced.engine->IterateTS(8);
        control.engine->IterateTS(8);
        for (unsigned int field = 0; field < 2; ++field)
        {
            const auto weight = std::polar(1e-14f, .2f * (batch + field));
            for (unsigned int run = 0; run < 2; ++run)
            {
                interfaces[run][field]->AccumulateFieldDFT(ids[run][field], {weight, -weight});
            }
        }
        const auto before = fenced.op->NativeObservationInfo();
        Reject([&] { fenced.op->FenceObservationEpoch(8 * (batch + 1) + 1); });
        for (unsigned int repeat = 0; repeat < 2; ++repeat)
        {
            Require(fenced.op->FenceObservationEpoch(8 * (batch + 1)) == workers,
                "native collective fence did not acknowledge every worker");
        }
        const auto after = fenced.op->NativeObservationInfo();
        Require(after.fullFieldReads == before.fullFieldReads && after.fieldDFTDownloadBytes == before.fieldDFTDownloadBytes &&
            after.selectedReadbacks == before.selectedReadbacks && after.hostObservationBytes == before.hostObservationBytes,
            "native fence downloaded or materialized fields/sums");
        for (unsigned int field = 0; field < 2; ++field)
        {
            for (unsigned int run = 0; run < 2; ++run)
            {
                std::vector<ArrayLib::ArrayNIJK<std::complex<float>> *> output{&values[run][0], &values[run][1]};
                Require(interfaces[run][field]->ReadFieldDFT(ids[run][field], output), "native fence DFT read failed");
            }
            for (unsigned int frequency = 0; frequency < 2; ++frequency)
            {
                for (std::size_t k = 0; k < values[0][frequency].size(); ++k)
                {
                    Require(values[0][frequency].data(k) == values[1][frequency].data(k), "fence changed native DFT arrays");
                }
            }
        }
    }
    const auto a = fenced.op->ObservationGrid().DeviceStateData(), b = control.op->ObservationGrid().DeviceStateData();
    Require(a.electric == b.electric && a.magnetic == b.magnetic, "native resume changed fields/material histories");
    Require(fenced.engine->GetNumberOfTimesteps() == 64 && control.engine->GetNumberOfTimesteps() == 64,
        "fence advanced or reset the native epoch");
    std::cout << "NATIVE_EPOCH_FENCE axis=" << axis << " loss=" << loss << " domains=" << domains
        << " workers=" << workers << " fences=17 exact_fields=1 exact_DFT=1 no_fence_download=1\n";
}

int main(int argc, char **argv)
{
    try
    {
        Require(argc <= 2, "unexpected native fence argument");
        const auto target = argc == 1 ? SubgridDeviceTarget{} : SubgridDeviceTarget::FromDomainFlag(argv[1]);
        for (unsigned int axis = 0; axis < 3; ++axis)
        {
            for (bool loss : {false, true})
            {
                for (bool domains : {false, true})
                {
                    Control(target, axis, loss, domains);
                }
            }
        }
        std::cout << "PASS: exact native epoch fences and uninterrupted fields/DFT after resume\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
