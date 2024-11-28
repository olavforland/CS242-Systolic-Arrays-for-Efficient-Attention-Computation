// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Design implementation internals
// See VattentionSystolicArray.h for the primary calling header

#include "VattentionSystolicArray__pch.h"
#include "VattentionSystolicArray_systolic_array__N4_M4.h"

VL_ATTR_COLD void VattentionSystolicArray_systolic_array__N4_M4___eval_initial__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray(VattentionSystolicArray_systolic_array__N4_M4* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+        VattentionSystolicArray_systolic_array__N4_M4___eval_initial__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelfRef.__PVT__partial_sum_wire[0U][0U] = 0.0;
    vlSelfRef.__PVT__partial_sum_wire[1U][0U] = 0.0;
    vlSelfRef.__PVT__partial_sum_wire[2U][0U] = 0.0;
    vlSelfRef.__PVT__partial_sum_wire[3U][0U] = 0.0;
}

VL_ATTR_COLD void VattentionSystolicArray_systolic_array__N4_M4___stl_sequent__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0(VattentionSystolicArray_systolic_array__N4_M4* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+        VattentionSystolicArray_systolic_array__N4_M4___stl_sequent__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelfRef.__PVT__data_wire[0U][0U] = vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[0U][1U] = vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[0U][2U] = vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[0U][3U] = vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[1U][0U] = vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[1U][1U] = vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[1U][2U] = vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[1U][3U] = vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[2U][0U] = vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[2U][1U] = vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[2U][2U] = vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[2U][3U] = vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[3U][0U] = vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[3U][1U] = vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[3U][2U] = vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__data_wire[3U][3U] = vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__data_out;
    vlSelfRef.__PVT__partial_sum_wire[0U][1U] = vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[0U][2U] = vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[0U][3U] = vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[0U][4U] = vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[1U][1U] = vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[1U][2U] = vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[1U][3U] = vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[1U][4U] = vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[2U][1U] = vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[2U][2U] = vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[2U][3U] = vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[2U][4U] = vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[3U][1U] = vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[3U][2U] = vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[3U][3U] = vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__partial_sum_wire[3U][4U] = vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out;
    vlSelfRef.__PVT__data_wire[4U][0U] = vlSelfRef.data_in
        [0U];
    vlSelfRef.__PVT__data_wire[4U][1U] = vlSelfRef.data_in
        [1U];
    vlSelfRef.__PVT__data_wire[4U][2U] = vlSelfRef.data_in
        [2U];
    vlSelfRef.__PVT__data_wire[4U][3U] = vlSelfRef.data_in
        [3U];
    vlSelfRef.result_out[0U] = vlSelfRef.__PVT__partial_sum_wire
        [0U][4U];
    vlSelfRef.result_out[1U] = vlSelfRef.__PVT__partial_sum_wire
        [1U][4U];
    vlSelfRef.result_out[2U] = vlSelfRef.__PVT__partial_sum_wire
        [2U][4U];
    vlSelfRef.result_out[3U] = vlSelfRef.__PVT__partial_sum_wire
        [3U][4U];
    vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[1U][0U] * vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[0U][0U]);
    vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[1U][1U] * vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[0U][1U]);
    vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[1U][2U] * vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[0U][2U]);
    vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[1U][3U] * vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[0U][3U]);
    vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[2U][0U] * vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[1U][0U]);
    vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[2U][1U] * vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[1U][1U]);
    vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[2U][2U] * vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[1U][2U]);
    vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[2U][3U] * vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[1U][3U]);
    vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[3U][0U] * vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[2U][0U]);
    vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[3U][1U] * vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[2U][1U]);
    vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[3U][2U] * vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[2U][2U]);
    vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[3U][3U] * vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[2U][3U]);
    vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[4U][0U] * vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[3U][0U]);
    vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[4U][1U] * vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[3U][1U]);
    vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[4U][2U] * vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[3U][2U]);
    vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result 
        = (vlSelfRef.__PVT__data_wire[4U][3U] * vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight);
    vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result 
        = (vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result 
           + vlSelfRef.__PVT__partial_sum_wire[3U][3U]);
}

VL_ATTR_COLD void VattentionSystolicArray_systolic_array__N4_M4___ctor_var_reset(VattentionSystolicArray_systolic_array__N4_M4* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+        VattentionSystolicArray_systolic_array__N4_M4___ctor_var_reset\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelf->clk = VL_RAND_RESET_I(1);
    vlSelf->reset = VL_RAND_RESET_I(1);
    vlSelf->weight_load_enable = VL_RAND_RESET_I(1);
    vlSelf->doProcess = VL_RAND_RESET_I(1);
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->weight_in[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        vlSelf->data_in[__Vi0] = 0;
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        vlSelf->result_out[__Vi0] = 0;
    }
    for (int __Vi0 = 0; __Vi0 < 5; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->__PVT__data_wire[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 5; ++__Vi1) {
            vlSelf->__PVT__partial_sum_wire[__Vi0][__Vi1] = 0;
        }
    }
    vlSelf->__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0;
    vlSelf->__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0;
    vlSelf->__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result = 0;
    vlSelf->__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result = 0;
}
