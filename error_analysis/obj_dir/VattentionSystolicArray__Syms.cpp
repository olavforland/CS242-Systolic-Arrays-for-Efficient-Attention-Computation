// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Symbol table implementation internals

#include "VattentionSystolicArray__pch.h"
#include "VattentionSystolicArray.h"
#include "VattentionSystolicArray___024root.h"
#include "VattentionSystolicArray_systolic_array__N4_M4.h"

// FUNCTIONS
VattentionSystolicArray__Syms::~VattentionSystolicArray__Syms()
{
}

VattentionSystolicArray__Syms::VattentionSystolicArray__Syms(VerilatedContext* contextp, const char* namep, VattentionSystolicArray* modelp)
    : VerilatedSyms{contextp}
    // Setup internal state of the Syms class
    , __Vm_modelp{modelp}
    // Setup module instances
    , TOP{this, namep}
    , TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray{this, Verilated::catName(namep, "attentionSystolicArray.K_mult_Q_systolicArray")}
    , TOP__attentionSystolicArray__DOT__V_mult_systolicArray{this, Verilated::catName(namep, "attentionSystolicArray.V_mult_systolicArray")}
{
        // Check resources
        Verilated::stackCheck(1247);
    // Configure time unit / time precision
    _vm_contextp__->timeunit(-9);
    _vm_contextp__->timeprecision(-12);
    // Setup each module's pointers to their submodules
    TOP.__PVT__attentionSystolicArray__DOT__K_mult_Q_systolicArray = &TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray;
    TOP.__PVT__attentionSystolicArray__DOT__V_mult_systolicArray = &TOP__attentionSystolicArray__DOT__V_mult_systolicArray;
    // Setup each module's pointer back to symbol table (for public functions)
    TOP.__Vconfigure(true);
    TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vconfigure(true);
    TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vconfigure(false);
}
