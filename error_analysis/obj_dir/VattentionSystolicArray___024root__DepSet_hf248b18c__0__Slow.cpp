// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Design implementation internals
// See VattentionSystolicArray.h for the primary calling header

#include "VattentionSystolicArray__pch.h"
#include "VattentionSystolicArray___024root.h"

VL_ATTR_COLD void VattentionSystolicArray___024root___eval_static(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_static\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
}

VL_ATTR_COLD void VattentionSystolicArray___024root___eval_initial__TOP(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_initial__TOP\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[0U][0U] = 1.0;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0U] = 1.0;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[1U][0U] = 1.0;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0U] = 1.0;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[2U][0U] = 1.0;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0U] = 1.0;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[3U][0U] = 1.0;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0U] = 1.0;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[4U] = 0.0;
}

VL_ATTR_COLD void VattentionSystolicArray___024root___eval_final(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_final\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
}

#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__stl(VattentionSystolicArray___024root* vlSelf);
#endif  // VL_DEBUG
VL_ATTR_COLD bool VattentionSystolicArray___024root___eval_phase__stl(VattentionSystolicArray___024root* vlSelf);

VL_ATTR_COLD void VattentionSystolicArray___024root___eval_settle(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_settle\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    IData/*31:0*/ __VstlIterCount;
    CData/*0:0*/ __VstlContinue;
    // Body
    __VstlIterCount = 0U;
    vlSelfRef.__VstlFirstIteration = 1U;
    __VstlContinue = 1U;
    while (__VstlContinue) {
        if (VL_UNLIKELY((0x64U < __VstlIterCount))) {
#ifdef VL_DEBUG
            VattentionSystolicArray___024root___dump_triggers__stl(vlSelf);
#endif
            VL_FATAL_MT("../rtl/../rtl/attentionSystolicArray.sv", 5, "", "Settle region did not converge.");
        }
        __VstlIterCount = ((IData)(1U) + __VstlIterCount);
        __VstlContinue = 0U;
        if (VattentionSystolicArray___024root___eval_phase__stl(vlSelf)) {
            __VstlContinue = 1U;
        }
        vlSelfRef.__VstlFirstIteration = 0U;
    }
}

#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__stl(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___dump_triggers__stl\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    if ((1U & (~ vlSelfRef.__VstlTriggered.any()))) {
        VL_DBG_MSGF("         No triggers active\n");
    }
    if ((1ULL & vlSelfRef.__VstlTriggered.word(0U))) {
        VL_DBG_MSGF("         'stl' region trigger index 0 is active: Internal 'stl' trigger - first iteration\n");
    }
}
#endif  // VL_DEBUG

VL_ATTR_COLD void VattentionSystolicArray___024root___eval_triggers__stl(VattentionSystolicArray___024root* vlSelf);
VL_ATTR_COLD void VattentionSystolicArray___024root___eval_stl(VattentionSystolicArray___024root* vlSelf);

VL_ATTR_COLD bool VattentionSystolicArray___024root___eval_phase__stl(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_phase__stl\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    CData/*0:0*/ __VstlExecute;
    // Body
    VattentionSystolicArray___024root___eval_triggers__stl(vlSelf);
    __VstlExecute = vlSelfRef.__VstlTriggered.any();
    if (__VstlExecute) {
        VattentionSystolicArray___024root___eval_stl(vlSelf);
    }
    return (__VstlExecute);
}

#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__ico(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___dump_triggers__ico\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    if ((1U & (~ vlSelfRef.__VicoTriggered.any()))) {
        VL_DBG_MSGF("         No triggers active\n");
    }
    if ((1ULL & vlSelfRef.__VicoTriggered.word(0U))) {
        VL_DBG_MSGF("         'ico' region trigger index 0 is active: Internal 'ico' trigger - first iteration\n");
    }
}
#endif  // VL_DEBUG

#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__act(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___dump_triggers__act\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    if ((1U & (~ vlSelfRef.__VactTriggered.any()))) {
        VL_DBG_MSGF("         No triggers active\n");
    }
    if ((1ULL & vlSelfRef.__VactTriggered.word(0U))) {
        VL_DBG_MSGF("         'act' region trigger index 0 is active: @(posedge clk)\n");
    }
}
#endif  // VL_DEBUG

#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__nba(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___dump_triggers__nba\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    if ((1U & (~ vlSelfRef.__VnbaTriggered.any()))) {
        VL_DBG_MSGF("         No triggers active\n");
    }
    if ((1ULL & vlSelfRef.__VnbaTriggered.word(0U))) {
        VL_DBG_MSGF("         'nba' region trigger index 0 is active: @(posedge clk)\n");
    }
}
#endif  // VL_DEBUG

VL_ATTR_COLD void VattentionSystolicArray___024root____Vm_traceActivitySetAll(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root____Vm_traceActivitySetAll\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelfRef.__Vm_traceActivity[0U] = 1U;
    vlSelfRef.__Vm_traceActivity[1U] = 1U;
    vlSelfRef.__Vm_traceActivity[2U] = 1U;
}

VL_ATTR_COLD void VattentionSystolicArray___024root___ctor_var_reset(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___ctor_var_reset\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelf->clk = VL_RAND_RESET_I(1);
    vlSelf->reset = VL_RAND_RESET_I(1);
    vlSelf->valid_input = VL_RAND_RESET_I(1);
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->K_matrix[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 33; ++__Vi0) {
        vlSelf->factorial_arr[__Vi0] = 0;
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->Q_matrix[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->V_matrix[__Vi0][__Vi1] = 0;
        }
    }
    vlSelf->softmax_temp = 0;
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->Q_mult_K[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->exp_Q_mult_K[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->exp_K_mult_Q_mult_V[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->attention[__Vi0][__Vi1] = 0;
        }
    }
    vlSelf->valid_result = VL_RAND_RESET_I(1);
    vlSelf->attentionSystolicArray__DOT__counter_q = 0;
    vlSelf->attentionSystolicArray__DOT__counter_d = 0;
    vlSelf->attentionSystolicArray__DOT__valid_result_q = VL_RAND_RESET_I(1);
    vlSelf->attentionSystolicArray__DOT__doProcess_d = VL_RAND_RESET_I(1);
    vlSelf->attentionSystolicArray__DOT__doProcess_q = VL_RAND_RESET_I(1);
    vlSelf->attentionSystolicArray__DOT__weight_load_enable_q = VL_RAND_RESET_I(1);
    vlSelf->attentionSystolicArray__DOT__weight_load_counter = VL_RAND_RESET_I(2);
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        vlSelf->attentionSystolicArray__DOT__data_in_array[__Vi0] = 0;
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->attentionSystolicArray__DOT__K_matrix_in[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        vlSelf->attentionSystolicArray__DOT__Q_mult_K_result[__Vi0] = 0;
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        vlSelf->attentionSystolicArray__DOT__exp_Q_mult_K_result[__Vi0] = 0;
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        vlSelf->attentionSystolicArray__DOT__Q_mult_K_result_normed[__Vi0] = 0;
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        vlSelf->attentionSystolicArray__DOT__data_in_v[__Vi0] = 0;
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 4; ++__Vi1) {
            vlSelf->attentionSystolicArray__DOT__weight_V_matrix[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        vlSelf->attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result[__Vi0] = 0;
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        vlSelf->attentionSystolicArray__DOT__attention_result[__Vi0] = 0;
    }
    vlSelf->attentionSystolicArray__DOT__unnamedblk1__DOT__j = 0;
    vlSelf->attentionSystolicArray__DOT__unnamedblk2__DOT__j = 0;
    vlSelf->attentionSystolicArray__DOT__unnamedblk9__DOT__i = 0;
    vlSelf->attentionSystolicArray__DOT__unnamedblk9__DOT__unnamedblk10__DOT__j = 0;
    vlSelf->attentionSystolicArray__DOT__unnamedblk11__DOT__i = 0;
    vlSelf->attentionSystolicArray__DOT__unnamedblk11__DOT__unnamedblk12__DOT__j = 0;
    vlSelf->attentionSystolicArray__DOT__unnamedblk13__DOT__i = 0;
    vlSelf->attentionSystolicArray__DOT__unnamedblk13__DOT__unnamedblk14__DOT__j = 0;
    vlSelf->attentionSystolicArray__DOT__unnamedblk15__DOT__i = 0;
    vlSelf->attentionSystolicArray__DOT__unnamedblk15__DOT__unnamedblk16__DOT__j = 0;
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 33; ++__Vi1) {
            vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 33; ++__Vi1) {
            vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 4; ++__Vi0) {
        for (int __Vi1 = 0; __Vi1 < 33; ++__Vi1) {
            vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[__Vi0][__Vi1] = 0;
        }
    }
    for (int __Vi0 = 0; __Vi0 < 5; ++__Vi0) {
        vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[__Vi0] = 0;
    }
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_sum_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_sum_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_sum_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__x_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_sum_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_out = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term = 0;
    vlSelf->attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum = 0;
    for (int __Vi0 = 0; __Vi0 < 5; ++__Vi0) {
        vlSelf->attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[__Vi0] = 0;
    }
    vlSelf->attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__dummy_pe__data_out = 0;
    vlSelf->attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__norm_out = 0;
    vlSelf->attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__x_normed_out = 0;
    vlSelf->attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__norm_out = 0;
    vlSelf->attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__x_normed_out = 0;
    vlSelf->attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__norm_out = 0;
    vlSelf->attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__x_normed_out = 0;
    vlSelf->attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__norm_out = 0;
    vlSelf->attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__x_normed_out = 0;
    vlSelf->__Vdly__attentionSystolicArray__DOT__weight_load_enable_q = VL_RAND_RESET_I(1);
    vlSelf->__Vtrigprevexpr___TOP__clk__0 = VL_RAND_RESET_I(1);
    for (int __Vi0 = 0; __Vi0 < 3; ++__Vi0) {
        vlSelf->__Vm_traceActivity[__Vi0] = 0;
    }
}
