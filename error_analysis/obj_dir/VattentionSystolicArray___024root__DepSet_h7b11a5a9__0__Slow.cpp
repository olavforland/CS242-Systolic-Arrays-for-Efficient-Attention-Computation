// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Design implementation internals
// See VattentionSystolicArray.h for the primary calling header

#include "VattentionSystolicArray__pch.h"
#include "VattentionSystolicArray__Syms.h"
#include "VattentionSystolicArray___024root.h"

VL_ATTR_COLD void VattentionSystolicArray___024root___eval_initial__TOP(VattentionSystolicArray___024root* vlSelf);
VL_ATTR_COLD void VattentionSystolicArray___024root____Vm_traceActivitySetAll(VattentionSystolicArray___024root* vlSelf);
VL_ATTR_COLD void VattentionSystolicArray_systolic_array__N4_M4___eval_initial__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray(VattentionSystolicArray_systolic_array__N4_M4* vlSelf);

VL_ATTR_COLD void VattentionSystolicArray___024root___eval_initial(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_initial\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    VattentionSystolicArray___024root___eval_initial__TOP(vlSelf);
    VattentionSystolicArray___024root____Vm_traceActivitySetAll(vlSelf);
    VattentionSystolicArray_systolic_array__N4_M4___eval_initial__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray((&vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray));
    VattentionSystolicArray_systolic_array__N4_M4___eval_initial__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray((&vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray));
    vlSelfRef.__Vtrigprevexpr___TOP__clk__0 = vlSelfRef.clk;
}

#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__stl(VattentionSystolicArray___024root* vlSelf);
#endif  // VL_DEBUG

VL_ATTR_COLD void VattentionSystolicArray___024root___eval_triggers__stl(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_triggers__stl\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelfRef.__VstlTriggered.set(0U, (IData)(vlSelfRef.__VstlFirstIteration));
#ifdef VL_DEBUG
    if (VL_UNLIKELY(vlSymsp->_vm_contextp__->debug())) {
        VattentionSystolicArray___024root___dump_triggers__stl(vlSelf);
    }
#endif
}

VL_ATTR_COLD void VattentionSystolicArray___024root___stl_sequent__TOP__0(VattentionSystolicArray___024root* vlSelf);
VL_ATTR_COLD void VattentionSystolicArray_systolic_array__N4_M4___stl_sequent__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0(VattentionSystolicArray_systolic_array__N4_M4* vlSelf);
VL_ATTR_COLD void VattentionSystolicArray___024root___stl_sequent__TOP__1(VattentionSystolicArray___024root* vlSelf);

VL_ATTR_COLD void VattentionSystolicArray___024root___eval_stl(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_stl\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    if ((1ULL & vlSelfRef.__VstlTriggered.word(0U))) {
        VattentionSystolicArray___024root___stl_sequent__TOP__0(vlSelf);
        VattentionSystolicArray___024root____Vm_traceActivitySetAll(vlSelf);
        VattentionSystolicArray_systolic_array__N4_M4___stl_sequent__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0((&vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray));
        VattentionSystolicArray_systolic_array__N4_M4___stl_sequent__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0((&vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray));
        VattentionSystolicArray___024root___stl_sequent__TOP__1(vlSelf);
    }
}

VL_ATTR_COLD void VattentionSystolicArray___024root___stl_sequent__TOP__0(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___stl_sequent__TOP__0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelfRef.valid_result = vlSelfRef.attentionSystolicArray__DOT__valid_result_q;
    vlSelfRef.attentionSystolicArray__DOT__counter_d 
        = ((IData)(vlSelfRef.attentionSystolicArray__DOT__doProcess_q)
            ? ((IData)(1U) + vlSelfRef.attentionSystolicArray__DOT__counter_q)
            : 0U);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[0U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_sum_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_sum_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_sum_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_sum_out;
    vlSelfRef.attentionSystolicArray__DOT__attention_result[0U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__x_normed_out;
    vlSelfRef.attentionSystolicArray__DOT__attention_result[1U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__x_normed_out;
    vlSelfRef.attentionSystolicArray__DOT__attention_result[2U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__x_normed_out;
    vlSelfRef.attentionSystolicArray__DOT__attention_result[3U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__x_normed_out;
    vlSelfRef.attentionSystolicArray__DOT__doProcess_d 
        = (((~ (IData)(vlSelfRef.attentionSystolicArray__DOT__weight_load_enable_q)) 
            & (IData)(vlSelfRef.valid_input)) | ((0x32U 
                                                  != vlSelfRef.attentionSystolicArray__DOT__counter_q) 
                                                 & (IData)(vlSelfRef.attentionSystolicArray__DOT__doProcess_q)));
    vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[4U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__dummy_pe__data_out;
    vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[0U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__norm_out;
    vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[1U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__norm_out;
    vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[2U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__norm_out;
    vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[3U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__norm_out;
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[0U][0U] 
        = vlSelfRef.K_matrix[3U][0U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[0U][1U] 
        = vlSelfRef.K_matrix[3U][1U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[0U][2U] 
        = vlSelfRef.K_matrix[3U][2U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[0U][3U] 
        = vlSelfRef.K_matrix[3U][3U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[1U][0U] 
        = vlSelfRef.K_matrix[2U][0U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[1U][1U] 
        = vlSelfRef.K_matrix[2U][1U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[1U][2U] 
        = vlSelfRef.K_matrix[2U][2U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[1U][3U] 
        = vlSelfRef.K_matrix[2U][3U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[2U][0U] 
        = vlSelfRef.K_matrix[1U][0U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[2U][1U] 
        = vlSelfRef.K_matrix[1U][1U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[2U][2U] 
        = vlSelfRef.K_matrix[1U][2U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[2U][3U] 
        = vlSelfRef.K_matrix[1U][3U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[3U][0U] 
        = vlSelfRef.K_matrix[0U][0U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[3U][1U] 
        = vlSelfRef.K_matrix[0U][1U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[3U][2U] 
        = vlSelfRef.K_matrix[0U][2U];
    vlSelfRef.attentionSystolicArray__DOT__K_matrix_in[3U][3U] 
        = vlSelfRef.K_matrix[0U][3U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[0U][0U] 
        = vlSelfRef.V_matrix[0U][3U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[0U][1U] 
        = vlSelfRef.V_matrix[1U][3U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[0U][2U] 
        = vlSelfRef.V_matrix[2U][3U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[0U][3U] 
        = vlSelfRef.V_matrix[3U][3U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[1U][0U] 
        = vlSelfRef.V_matrix[0U][2U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[1U][1U] 
        = vlSelfRef.V_matrix[1U][2U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[1U][2U] 
        = vlSelfRef.V_matrix[2U][2U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[1U][3U] 
        = vlSelfRef.V_matrix[3U][2U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[2U][0U] 
        = vlSelfRef.V_matrix[0U][1U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[2U][1U] 
        = vlSelfRef.V_matrix[1U][1U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[2U][2U] 
        = vlSelfRef.V_matrix[2U][1U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[2U][3U] 
        = vlSelfRef.V_matrix[3U][1U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[3U][0U] 
        = vlSelfRef.V_matrix[0U][0U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[3U][1U] 
        = vlSelfRef.V_matrix[1U][0U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[3U][2U] 
        = vlSelfRef.V_matrix[2U][0U];
    vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix[3U][3U] 
        = vlSelfRef.V_matrix[3U][0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.data_in[0U] 
        = vlSelfRef.attentionSystolicArray__DOT__data_in_array
        [0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.data_in[1U] 
        = vlSelfRef.attentionSystolicArray__DOT__data_in_array
        [1U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.data_in[2U] 
        = vlSelfRef.attentionSystolicArray__DOT__data_in_array
        [2U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.data_in[3U] 
        = vlSelfRef.attentionSystolicArray__DOT__data_in_array
        [3U];
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__x_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x20U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x20U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x20U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][4U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][5U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][6U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][7U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][8U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][9U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0xaU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0xbU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0xcU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0xdU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0xeU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0xfU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x10U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x11U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x12U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x13U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x14U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x15U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x16U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x17U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x18U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x19U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x1aU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x1bU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x1cU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x1dU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x1eU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x1fU] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x20U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_out;
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[0U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [0U][0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[0U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [0U][1U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[0U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [0U][2U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[0U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [0U][3U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[1U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [1U][0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[1U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [1U][1U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[1U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [1U][2U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[1U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [1U][3U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[2U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [2U][0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[2U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [2U][1U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[2U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [2U][2U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[2U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [2U][3U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[3U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [3U][0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[3U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [3U][1U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[3U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [3U][2U];
    vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in[3U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
        [3U][3U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[0U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [0U][0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[0U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [0U][1U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[0U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [0U][2U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[0U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [0U][3U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[1U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [1U][0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[1U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [1U][1U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[1U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [1U][2U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[1U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [1U][3U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[2U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [2U][0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[2U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [2U][1U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[2U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [2U][2U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[2U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [2U][3U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[3U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [3U][0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[3U][1U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [3U][1U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[3U][2U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [3U][2U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in[3U][3U] 
        = vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
        [3U][3U];
    vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result[0U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
        [0U][0x20U];
    vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result[1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
        [1U][0x20U];
    vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result[2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
        [2U][0x20U];
    vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result[3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
        [3U][0x20U];
    vlSelfRef.attentionSystolicArray__DOT__data_in_v[0U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
        [3U];
    vlSelfRef.attentionSystolicArray__DOT__data_in_v[1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
        [2U];
    vlSelfRef.attentionSystolicArray__DOT__data_in_v[2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
        [1U];
    vlSelfRef.attentionSystolicArray__DOT__data_in_v[3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
        [0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.data_in[0U] 
        = vlSelfRef.attentionSystolicArray__DOT__data_in_v
        [0U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.data_in[1U] 
        = vlSelfRef.attentionSystolicArray__DOT__data_in_v
        [1U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.data_in[2U] 
        = vlSelfRef.attentionSystolicArray__DOT__data_in_v
        [2U];
    vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.data_in[3U] 
        = vlSelfRef.attentionSystolicArray__DOT__data_in_v
        [3U];
}

VL_ATTR_COLD void VattentionSystolicArray___024root___stl_sequent__TOP__1(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___stl_sequent__TOP__1\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result[0U] 
        = vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.result_out
        [0U];
    vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result[1U] 
        = vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.result_out
        [1U];
    vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result[2U] 
        = vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.result_out
        [2U];
    vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result[3U] 
        = vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.result_out
        [3U];
    vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result[0U] 
        = vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.result_out
        [0U];
    vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result[1U] 
        = vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.result_out
        [1U];
    vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result[2U] 
        = vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.result_out
        [2U];
    vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result[3U] 
        = vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.result_out
        [3U];
    vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed[0U] 
        = (vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
           [0U] * vlSelfRef.softmax_temp);
    vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed[1U] 
        = (vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
           [1U] * vlSelfRef.softmax_temp);
    vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed[2U] 
        = (vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
           [2U] * vlSelfRef.softmax_temp);
    vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed[3U] 
        = (vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
           [3U] * vlSelfRef.softmax_temp);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[0U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed
        [0U];
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[1U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed
        [1U];
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[2U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed
        [2U];
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[3U][0U] 
        = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed
        [3U];
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[1U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][1U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][1U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[2U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][1U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][2U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][2U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[3U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][2U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][3U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][3U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[4U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][3U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][4U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][4U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[5U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][4U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][5U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][5U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[6U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][5U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][6U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][6U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[7U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][6U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][7U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][7U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[8U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][7U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][8U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][8U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[9U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][8U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][9U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][9U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xaU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][9U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0xaU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0xaU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xbU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0xaU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0xbU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0xbU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xcU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0xbU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0xcU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0xcU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xdU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0xcU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0xdU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0xdU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xeU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0xdU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0xeU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0xeU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xfU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0xeU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0xfU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0xfU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x10U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0xfU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x10U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x10U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x11U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x10U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x11U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x11U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x12U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x11U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x12U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x12U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x13U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x12U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x13U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x13U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x14U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x13U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x14U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x14U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x15U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x14U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x15U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x15U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x16U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x15U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x16U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x16U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x17U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x16U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x17U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x17U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x18U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x17U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x18U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x18U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x19U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x18U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x19U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x19U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1aU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x19U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x1aU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x1aU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1bU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x1aU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x1bU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x1bU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1cU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x1bU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x1cU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x1cU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1dU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x1cU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x1dU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x1dU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1eU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x1dU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [0U][0x1eU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [0U][0x1eU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1fU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [0U][0x1eU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[1U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][1U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][1U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[2U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][1U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][2U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][2U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[3U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][2U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][3U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][3U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[4U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][3U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][4U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][4U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[5U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][4U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][5U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][5U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[6U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][5U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][6U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][6U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[7U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][6U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][7U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][7U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[8U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][7U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][8U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][8U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[9U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][8U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][9U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][9U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xaU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][9U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0xaU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0xaU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xbU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0xaU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0xbU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0xbU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xcU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0xbU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0xcU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0xcU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xdU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0xcU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0xdU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0xdU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xeU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0xdU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0xeU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0xeU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xfU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0xeU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0xfU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0xfU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x10U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0xfU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x10U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x10U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x11U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x10U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x11U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x11U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x12U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x11U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x12U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x12U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x13U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x12U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x13U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x13U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x14U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x13U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x14U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x14U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x15U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x14U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x15U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x15U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x16U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x15U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x16U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x16U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x17U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x16U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x17U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x17U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x18U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x17U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x18U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x18U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x19U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x18U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x19U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x19U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1aU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x19U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x1aU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x1aU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1bU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x1aU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x1bU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x1bU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1cU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x1bU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x1cU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x1cU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1dU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x1cU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x1dU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x1dU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1eU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x1dU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [1U][0x1eU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [1U][0x1eU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1fU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [1U][0x1eU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[1U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][1U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][1U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[2U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][1U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][2U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][2U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[3U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][2U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][3U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][3U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[4U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][3U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][4U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][4U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[5U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][4U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][5U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][5U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[6U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][5U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][6U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][6U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[7U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][6U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][7U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][7U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[8U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][7U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][8U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][8U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[9U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][8U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][9U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][9U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xaU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][9U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0xaU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0xaU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xbU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0xaU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0xbU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0xbU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xcU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0xbU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0xcU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0xcU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xdU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0xcU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0xdU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0xdU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xeU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0xdU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0xeU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0xeU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xfU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0xeU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0xfU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0xfU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x10U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0xfU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x10U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x10U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x11U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x10U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x11U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x11U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x12U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x11U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x12U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x12U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x13U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x12U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x13U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x13U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x14U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x13U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x14U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x14U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x15U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x14U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x15U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x15U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x16U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x15U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x16U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x16U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x17U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x16U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x17U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x17U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x18U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x17U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x18U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x18U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x19U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x18U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x19U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x19U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1aU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x19U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x1aU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x1aU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1bU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x1aU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x1bU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x1bU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1cU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x1bU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x1cU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x1cU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1dU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x1cU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x1dU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x1dU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1eU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x1dU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [2U][0x1eU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [2U][0x1eU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1fU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [2U][0x1eU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[1U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][1U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][1U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[2U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][1U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][2U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][2U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[3U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][2U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][3U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][3U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[4U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][3U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][4U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][4U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[5U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][4U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][5U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][5U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[6U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][5U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][6U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][6U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[7U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][6U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][7U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][7U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[8U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][7U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][8U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][8U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[9U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][8U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][9U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][9U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xaU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][9U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0xaU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0xaU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xbU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0xaU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0xbU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0xbU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xcU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0xbU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0xcU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0xcU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xdU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0xcU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0xdU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0xdU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xeU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0xdU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0xeU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0xeU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0xfU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0xeU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0xfU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0xfU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x10U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0xfU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x10U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x10U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x11U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x10U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x11U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x11U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x12U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x11U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x12U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x12U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x13U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x12U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x13U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x13U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x14U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x13U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x14U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x14U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x15U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x14U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x15U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x15U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x16U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x15U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x16U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x16U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x17U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x16U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x17U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x17U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x18U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x17U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x18U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x18U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x19U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x18U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x19U] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x19U]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1aU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x19U] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x1aU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x1aU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1bU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x1aU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x1bU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x1bU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1cU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x1bU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x1cU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x1cU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1dU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x1cU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x1dU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x1dU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1eU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x1dU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
           [3U][0x1eU] * vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
           [3U][0x1eU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term 
           / vlSelfRef.factorial_arr[0x1fU]);
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum 
        = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
           [3U][0x1eU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term);
}
