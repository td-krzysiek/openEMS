// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef OPENEMS_DISTRIBUTED_OPERATOR_SLICE_H
#define OPENEMS_DISTRIBUTED_OPERATOR_SLICE_H
#include "cpu_executor.h"
class Operator;

namespace openems
{
namespace distributed
{
// Integration boundary: preserve globally averaged, fully built coefficients.
// Dynamic extensions are rejected, never dropped or replaced with PEC.
LocalOperatorData SliceOperator(const Operator&, const Domain&);
}
}
#endif
