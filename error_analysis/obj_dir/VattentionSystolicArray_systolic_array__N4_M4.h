// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Design internal header
// See VattentionSystolicArray.h for the primary calling header

#ifndef VERILATED_VATTENTIONSYSTOLICARRAY_SYSTOLIC_ARRAY__N4_M4_H_
#define VERILATED_VATTENTIONSYSTOLICARRAY_SYSTOLIC_ARRAY__N4_M4_H_  // guard

#include "verilated.h"


class VattentionSystolicArray__Syms;

class alignas(VL_CACHE_LINE_BYTES) VattentionSystolicArray_systolic_array__N4_M4 final : public VerilatedModule {
  public:

    // DESIGN SPECIFIC STATE
    VL_IN8(clk,0,0);
    VL_IN8(reset,0,0);
    VL_IN8(weight_load_enable,0,0);
    VL_IN8(doProcess,0,0);
    double __Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__data_out;
    double __Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out;
    double __Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__data_out;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result;
    double __PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
    VlUnpacked<VlUnpacked<double, 4>, 4> weight_in;
    VlUnpacked<double, 4> data_in;
    VlUnpacked<double, 4> result_out;
    VlUnpacked<VlUnpacked<double, 4>, 5> __PVT__data_wire;
    VlUnpacked<VlUnpacked<double, 5>, 4> __PVT__partial_sum_wire;

    // INTERNAL VARIABLES
    VattentionSystolicArray__Syms* const vlSymsp;

    // CONSTRUCTORS
    VattentionSystolicArray_systolic_array__N4_M4(VattentionSystolicArray__Syms* symsp, const char* v__name);
    ~VattentionSystolicArray_systolic_array__N4_M4();
    VL_UNCOPYABLE(VattentionSystolicArray_systolic_array__N4_M4);

    // INTERNAL METHODS
    void __Vconfigure(bool first);
};


#endif  // guard
