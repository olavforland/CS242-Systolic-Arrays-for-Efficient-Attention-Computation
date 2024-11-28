// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Design implementation internals
// See VattentionSystolicArray.h for the primary calling header

#include "VattentionSystolicArray__pch.h"
#include "VattentionSystolicArray__Syms.h"
#include "VattentionSystolicArray_systolic_array__N4_M4.h"

VL_INLINE_OPT void VattentionSystolicArray_systolic_array__N4_M4___nba_sequent__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0(VattentionSystolicArray_systolic_array__N4_M4* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+        VattentionSystolicArray_systolic_array__N4_M4___nba_sequent__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    if (vlSymsp->TOP.reset) {
        vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0.0;
    } else {
        if (vlSymsp->TOP.attentionSystolicArray__DOT__weight_load_enable_q) {
            vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[0U][0U];
            vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[0U][1U];
            vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[0U][2U];
            vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[0U][3U];
            vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[1U][0U];
            vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[1U][1U];
            vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[1U][2U];
            vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[1U][3U];
            vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[2U][0U];
            vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[2U][1U];
            vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[2U][2U];
            vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[2U][3U];
            vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[3U][0U];
            vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[3U][1U];
            vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[3U][2U];
            vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[3U][3U];
        }
        if (vlSymsp->TOP.attentionSystolicArray__DOT__doProcess_q) {
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[1U][0U];
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[1U][1U];
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[1U][2U];
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[1U][3U];
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[2U][0U];
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[2U][1U];
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[2U][2U];
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[2U][3U];
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[3U][0U];
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[3U][1U];
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[3U][2U];
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[3U][3U];
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[4U][0U];
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[4U][1U];
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[4U][2U];
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[4U][3U];
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
        }
    }
    vlSelfRef.__PVT__data_wire[4U][0U] = vlSelfRef.data_in
        [0U];
    vlSelfRef.__PVT__data_wire[4U][1U] = vlSelfRef.data_in
        [1U];
    vlSelfRef.__PVT__data_wire[4U][2U] = vlSelfRef.data_in
        [2U];
    vlSelfRef.__PVT__data_wire[4U][3U] = vlSelfRef.data_in
        [3U];
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

VL_INLINE_OPT void VattentionSystolicArray_systolic_array__N4_M4___nba_sequent__TOP__attentionSystolicArray__DOT__V_mult_systolicArray__0(VattentionSystolicArray_systolic_array__N4_M4* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+        VattentionSystolicArray_systolic_array__N4_M4___nba_sequent__TOP__attentionSystolicArray__DOT__V_mult_systolicArray__0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    if (vlSymsp->TOP.reset) {
        vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__data_out = 0.0;
        vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__data_out = 0.0;
    } else {
        if (vlSymsp->TOP.attentionSystolicArray__DOT__weight_load_enable_q) {
            vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[0U][0U];
            vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[0U][1U];
            vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[0U][2U];
            vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[0U][3U];
            vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[1U][0U];
            vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[1U][1U];
            vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[1U][2U];
            vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[1U][3U];
            vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[2U][0U];
            vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[2U][1U];
            vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[2U][2U];
            vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[2U][3U];
            vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[3U][0U];
            vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[3U][1U];
            vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[3U][2U];
            vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight 
                = vlSelfRef.weight_in[3U][3U];
        }
        if (vlSymsp->TOP.attentionSystolicArray__DOT__doProcess_q) {
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out 
                = vlSelfRef.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result;
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[1U][0U];
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[1U][1U];
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[1U][2U];
            vlSelfRef.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[1U][3U];
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[2U][0U];
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[2U][1U];
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[2U][2U];
            vlSelfRef.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[2U][3U];
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[3U][0U];
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[3U][1U];
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[3U][2U];
            vlSelfRef.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[3U][3U];
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[4U][0U];
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[4U][1U];
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[4U][2U];
            vlSelfRef.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__data_out 
                = vlSelfRef.__PVT__data_wire[4U][3U];
        }
    }
    vlSelfRef.__PVT__data_wire[4U][0U] = vlSelfRef.data_in
        [0U];
    vlSelfRef.__PVT__data_wire[4U][1U] = vlSelfRef.data_in
        [1U];
    vlSelfRef.__PVT__data_wire[4U][2U] = vlSelfRef.data_in
        [2U];
    vlSelfRef.__PVT__data_wire[4U][3U] = vlSelfRef.data_in
        [3U];
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
