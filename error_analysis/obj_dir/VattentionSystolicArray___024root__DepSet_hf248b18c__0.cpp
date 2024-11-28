// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Design implementation internals
// See VattentionSystolicArray.h for the primary calling header

#include "VattentionSystolicArray__pch.h"
#include "VattentionSystolicArray___024root.h"

void VattentionSystolicArray___024root___ico_sequent__TOP__0(VattentionSystolicArray___024root* vlSelf);

void VattentionSystolicArray___024root___eval_ico(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_ico\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    if ((1ULL & vlSelfRef.__VicoTriggered.word(0U))) {
        VattentionSystolicArray___024root___ico_sequent__TOP__0(vlSelf);
        vlSelfRef.__Vm_traceActivity[1U] = 1U;
    }
}

void VattentionSystolicArray___024root___eval_triggers__ico(VattentionSystolicArray___024root* vlSelf);

bool VattentionSystolicArray___024root___eval_phase__ico(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_phase__ico\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    CData/*0:0*/ __VicoExecute;
    // Body
    VattentionSystolicArray___024root___eval_triggers__ico(vlSelf);
    __VicoExecute = vlSelfRef.__VicoTriggered.any();
    if (__VicoExecute) {
        VattentionSystolicArray___024root___eval_ico(vlSelf);
    }
    return (__VicoExecute);
}

void VattentionSystolicArray___024root___eval_act(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_act\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
}

void VattentionSystolicArray___024root___eval_triggers__act(VattentionSystolicArray___024root* vlSelf);

bool VattentionSystolicArray___024root___eval_phase__act(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_phase__act\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    VlTriggerVec<1> __VpreTriggered;
    CData/*0:0*/ __VactExecute;
    // Body
    VattentionSystolicArray___024root___eval_triggers__act(vlSelf);
    __VactExecute = vlSelfRef.__VactTriggered.any();
    if (__VactExecute) {
        __VpreTriggered.andNot(vlSelfRef.__VactTriggered, vlSelfRef.__VnbaTriggered);
        vlSelfRef.__VnbaTriggered.thisOr(vlSelfRef.__VactTriggered);
        VattentionSystolicArray___024root___eval_act(vlSelf);
    }
    return (__VactExecute);
}

void VattentionSystolicArray___024root___eval_nba(VattentionSystolicArray___024root* vlSelf);

bool VattentionSystolicArray___024root___eval_phase__nba(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_phase__nba\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    CData/*0:0*/ __VnbaExecute;
    // Body
    __VnbaExecute = vlSelfRef.__VnbaTriggered.any();
    if (__VnbaExecute) {
        VattentionSystolicArray___024root___eval_nba(vlSelf);
        vlSelfRef.__VnbaTriggered.clear();
    }
    return (__VnbaExecute);
}

#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__ico(VattentionSystolicArray___024root* vlSelf);
#endif  // VL_DEBUG
#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__nba(VattentionSystolicArray___024root* vlSelf);
#endif  // VL_DEBUG
#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__act(VattentionSystolicArray___024root* vlSelf);
#endif  // VL_DEBUG

void VattentionSystolicArray___024root___eval(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    IData/*31:0*/ __VicoIterCount;
    CData/*0:0*/ __VicoContinue;
    IData/*31:0*/ __VnbaIterCount;
    CData/*0:0*/ __VnbaContinue;
    // Body
    __VicoIterCount = 0U;
    vlSelfRef.__VicoFirstIteration = 1U;
    __VicoContinue = 1U;
    while (__VicoContinue) {
        if (VL_UNLIKELY((0x64U < __VicoIterCount))) {
#ifdef VL_DEBUG
            VattentionSystolicArray___024root___dump_triggers__ico(vlSelf);
#endif
            VL_FATAL_MT("../rtl/../rtl/attentionSystolicArray.sv", 5, "", "Input combinational region did not converge.");
        }
        __VicoIterCount = ((IData)(1U) + __VicoIterCount);
        __VicoContinue = 0U;
        if (VattentionSystolicArray___024root___eval_phase__ico(vlSelf)) {
            __VicoContinue = 1U;
        }
        vlSelfRef.__VicoFirstIteration = 0U;
    }
    __VnbaIterCount = 0U;
    __VnbaContinue = 1U;
    while (__VnbaContinue) {
        if (VL_UNLIKELY((0x64U < __VnbaIterCount))) {
#ifdef VL_DEBUG
            VattentionSystolicArray___024root___dump_triggers__nba(vlSelf);
#endif
            VL_FATAL_MT("../rtl/../rtl/attentionSystolicArray.sv", 5, "", "NBA region did not converge.");
        }
        __VnbaIterCount = ((IData)(1U) + __VnbaIterCount);
        __VnbaContinue = 0U;
        vlSelfRef.__VactIterCount = 0U;
        vlSelfRef.__VactContinue = 1U;
        while (vlSelfRef.__VactContinue) {
            if (VL_UNLIKELY((0x64U < vlSelfRef.__VactIterCount))) {
#ifdef VL_DEBUG
                VattentionSystolicArray___024root___dump_triggers__act(vlSelf);
#endif
                VL_FATAL_MT("../rtl/../rtl/attentionSystolicArray.sv", 5, "", "Active region did not converge.");
            }
            vlSelfRef.__VactIterCount = ((IData)(1U) 
                                         + vlSelfRef.__VactIterCount);
            vlSelfRef.__VactContinue = 0U;
            if (VattentionSystolicArray___024root___eval_phase__act(vlSelf)) {
                vlSelfRef.__VactContinue = 1U;
            }
        }
        if (VattentionSystolicArray___024root___eval_phase__nba(vlSelf)) {
            __VnbaContinue = 1U;
        }
    }
}

#ifdef VL_DEBUG
void VattentionSystolicArray___024root___eval_debug_assertions(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_debug_assertions\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    if (VL_UNLIKELY((vlSelfRef.clk & 0xfeU))) {
        Verilated::overWidthError("clk");}
    if (VL_UNLIKELY((vlSelfRef.reset & 0xfeU))) {
        Verilated::overWidthError("reset");}
    if (VL_UNLIKELY((vlSelfRef.valid_input & 0xfeU))) {
        Verilated::overWidthError("valid_input");}
}
#endif  // VL_DEBUG
