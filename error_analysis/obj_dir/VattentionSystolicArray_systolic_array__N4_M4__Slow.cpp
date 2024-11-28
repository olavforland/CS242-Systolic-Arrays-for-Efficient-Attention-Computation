// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Design implementation internals
// See VattentionSystolicArray.h for the primary calling header

#include "VattentionSystolicArray__pch.h"
#include "VattentionSystolicArray__Syms.h"
#include "VattentionSystolicArray_systolic_array__N4_M4.h"

void VattentionSystolicArray_systolic_array__N4_M4___ctor_var_reset(VattentionSystolicArray_systolic_array__N4_M4* vlSelf);

VattentionSystolicArray_systolic_array__N4_M4::VattentionSystolicArray_systolic_array__N4_M4(VattentionSystolicArray__Syms* symsp, const char* v__name)
    : VerilatedModule{v__name}
    , vlSymsp{symsp}
 {
    // Reset structure values
    VattentionSystolicArray_systolic_array__N4_M4___ctor_var_reset(this);
}

void VattentionSystolicArray_systolic_array__N4_M4::__Vconfigure(bool first) {
    (void)first;  // Prevent unused variable warning
}

VattentionSystolicArray_systolic_array__N4_M4::~VattentionSystolicArray_systolic_array__N4_M4() {
}
