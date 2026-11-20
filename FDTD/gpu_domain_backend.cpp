// SPDX-License-Identifier: GPL-3.0-or-later
#include "gpu_backend.h"
#include "gpu_cell_model.h"
#include "gpu_domain_layout.h"
#include "operator.h"

bool GPU_Backend::InitDomain(const Operator* op, const GPU_CellModel& model,
    unsigned int begin, unsigned int end)
{
    m_Domain.reset(new GPU_DomainLayout());
    m_Domain->enabled = true;
    m_Domain->begin = begin;
    m_Domain->end = end;
    std::copy(model.numLines, model.numLines + 3, m_Domain->global);
    GPU_CellModel local(model, begin, end);
    return Init(op, local);
}

Engine_Box GPU_Backend::OwnedBox(bool currents) const
{
    if (!m_Domain)
        throw std::logic_error("owned updates require an initialized domain");
    const auto& d = *m_Domain;
    Engine_Box box{{1, 0, 0}, {d.LocalLines()-1, d.global[1], d.global[2]}};
    if (currents)
    {
        --box.stop[1];
        --box.stop[2];
        if (d.end == d.global[0])
            --box.stop[0];
    }
    return box;
}

void GPU_Backend::UpdateCurrentRegion(bool interior)
{
    const Engine_Box all = OwnedBox(true);
    if (all.Empty()) return;
    Engine_Box middle = all;
    // Only the first/last owned planes depend on an incoming voltage halo.
    if (m_Domain->begin > 0)
        ++middle.start[0];
    if (m_Domain->end < m_Domain->global[0] && middle.stop[0] > middle.start[0])
        --middle.stop[0];
    if (interior)
    {
        if (!middle.Empty()) UpdateCurrents(middle);
        return;
    }
    Engine_Box lower = all, upper = all;
    lower.stop[0] = middle.start[0];
    upper.start[0] = middle.stop[0];
    if (!lower.Empty()) UpdateCurrents(lower);
    if (!upper.Empty()) UpdateCurrents(upper);
}
