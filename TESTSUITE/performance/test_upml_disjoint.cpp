// SPDX-License-Identifier: GPL-3.0-or-later
// Regression for the barrier-group safety predicate (experimental patch).
#ifdef NDEBUG
#error "Compile this regression with assertions enabled"
#endif
#include "FDTD/extensions/engine_ext_upml.h"
#include "FDTD/extensions/operator_ext_upml.h"
#include <array>
#include <cassert>
#include <limits>
#include <iostream>

class TestOperator : public Operator_Ext_UPML
{
public:
    TestOperator(std::array<unsigned int, 3> start, std::array<unsigned int, 3> stop)
        : Operator_Ext_UPML(NULL) {SetRange(start.data(), stop.data());}
};

class CustomUPML : public Engine_Ext_UPML
{
public:
    explicit CustomUPML(Operator_Ext_UPML* op) : Engine_Ext_UPML(op) {}
};

int main()
{
    TestOperator oa({0,0,0}, {2,2,2});
    Engine_Ext_UPML a(&oa);
    assert(!a.IsDisjoint(a));
    for (unsigned int axis=0; axis<3; ++axis)
    {
        std::array<unsigned int,3> start = {0,0,0}, stop = {2,2,2};
        start[axis] = 3; stop[axis] = 5;
        TestOperator ob(start, stop);
        Engine_Ext_UPML b(&ob);
        assert(a.IsDisjoint(b) && b.IsDisjoint(a));
        CustomUPML custom(&ob);
        assert(!a.IsDisjoint(custom) && !custom.IsDisjoint(a));

        start[axis] = 2; stop[axis] = 4; // shared plane must keep the barrier
        TestOperator oc(start, stop);
        Engine_Ext_UPML c(&oc);
        assert(!a.IsDisjoint(c) && !c.IsDisjoint(a));
    }
    TestOperator inside({1,1,1}, {1,1,1});
    Engine_Ext_UPML contained(&inside);
    assert(!a.IsDisjoint(contained) && !contained.IsDisjoint(a));
    const unsigned int maximum = std::numeric_limits<unsigned int>::max();
    TestOperator high({maximum-2,0,0}, {maximum,2,2});
    Engine_Ext_UPML at_limit(&high);
    assert(!at_limit.IsDisjoint(at_limit));
    assert(at_limit.IsDisjoint(a) && a.IsDisjoint(at_limit));
    std::cout << "UPML disjoint-region safety checks passed\n";
}
