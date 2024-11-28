// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Symbol table internal header
//
// Internal details; most calling programs do not need this header,
// unless using verilator public meta comments.

#ifndef VERILATED_VATTENTIONSYSTOLICARRAY__SYMS_H_
#define VERILATED_VATTENTIONSYSTOLICARRAY__SYMS_H_  // guard

#include "verilated.h"

// INCLUDE MODEL CLASS

#include "VattentionSystolicArray.h"

// INCLUDE MODULE CLASSES
#include "VattentionSystolicArray___024root.h"
#include "VattentionSystolicArray_systolic_array__N4_M4.h"

// SYMS CLASS (contains all model state)
class alignas(VL_CACHE_LINE_BYTES)VattentionSystolicArray__Syms final : public VerilatedSyms {
  public:
    // INTERNAL STATE
    VattentionSystolicArray* const __Vm_modelp;
    bool __Vm_activity = false;  ///< Used by trace routines to determine change occurred
    uint32_t __Vm_baseCode = 0;  ///< Used by trace routines when tracing multiple models
    VlDeleter __Vm_deleter;
    bool __Vm_didInit = false;

    // MODULE INSTANCE STATE
    VattentionSystolicArray___024root TOP;
    VattentionSystolicArray_systolic_array__N4_M4 TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray;
    VattentionSystolicArray_systolic_array__N4_M4 TOP__attentionSystolicArray__DOT__V_mult_systolicArray;

    // CONSTRUCTORS
    VattentionSystolicArray__Syms(VerilatedContext* contextp, const char* namep, VattentionSystolicArray* modelp);
    ~VattentionSystolicArray__Syms();

    // METHODS
    const char* name() { return TOP.name(); }
};

#endif  // guard
