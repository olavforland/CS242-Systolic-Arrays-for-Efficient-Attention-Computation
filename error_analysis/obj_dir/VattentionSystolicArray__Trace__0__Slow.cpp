// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Tracing implementation internals
#include "verilated_vcd_c.h"
#include "VattentionSystolicArray__Syms.h"


VL_ATTR_COLD void VattentionSystolicArray___024root__trace_init_sub__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0(VattentionSystolicArray___024root* vlSelf, VerilatedVcd* tracep);
VL_ATTR_COLD void VattentionSystolicArray___024root__trace_init_sub__TOP__attentionSystolicArray__DOT__V_mult_systolicArray__0(VattentionSystolicArray___024root* vlSelf, VerilatedVcd* tracep);

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_init_sub__TOP__0(VattentionSystolicArray___024root* vlSelf, VerilatedVcd* tracep) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root__trace_init_sub__TOP__0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    const int c = vlSymsp->__Vm_baseCode;
    // Body
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3245,0,"valid_input",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->pushPrefix("K_matrix", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3246,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3248,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3250,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3252,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3254,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3256,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3258,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3260,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3262,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3264,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3266,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3268,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3270,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3272,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3274,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3276,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("Q_matrix", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3278,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3280,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3282,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3284,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3286,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3288,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3290,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3292,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3294,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3296,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3298,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3300,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3302,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3304,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3306,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3308,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("V_matrix", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3310,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3312,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3314,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3316,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3318,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3320,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3322,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3324,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3326,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3328,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3330,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3332,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3334,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3336,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3338,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3340,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->declDouble(c+3342,0,"softmax_temp",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->pushPrefix("Q_mult_K", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3344,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3346,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3348,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3350,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3352,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3354,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3356,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3358,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3360,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3362,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3364,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3366,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3368,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3370,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3372,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3374,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("exp_Q_mult_K", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3376,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3378,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3380,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3382,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3384,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3386,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3388,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3390,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3392,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3394,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3396,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3398,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3400,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3402,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3404,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3406,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("exp_K_mult_Q_mult_V", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3408,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3410,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3412,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3414,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3416,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3418,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3420,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3422,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3424,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3426,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3428,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3430,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3432,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3434,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3436,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3438,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("attention", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3440,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3442,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3444,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3446,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3448,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3450,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3452,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3454,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3456,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3458,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3460,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3462,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3464,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3466,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3468,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3470,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->declBit(c+3472,0,"valid_result",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->pushPrefix("attentionSystolicArray", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+3535,0,"N",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBus(c+3535,0,"d",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBus(c+3536,0,"K",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3245,0,"valid_input",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->pushPrefix("K_matrix", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3246,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3248,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3250,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3252,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3254,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3256,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3258,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3260,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3262,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3264,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3266,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3268,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3270,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3272,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3274,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3276,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("Q_matrix", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3278,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3280,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3282,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3284,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3286,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3288,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3290,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3292,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3294,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3296,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3298,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3300,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3302,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3304,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3306,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3308,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("V_matrix", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3310,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3312,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3314,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3316,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3318,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3320,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3322,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3324,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3326,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3328,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3330,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3332,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3334,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3336,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3338,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3340,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->declDouble(c+3342,0,"softmax_temp",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->pushPrefix("Q_mult_K", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3344,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3346,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3348,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3350,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3352,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3354,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3356,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3358,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3360,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3362,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3364,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3366,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3368,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3370,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3372,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3374,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("exp_Q_mult_K", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3376,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3378,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3380,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3382,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3384,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3386,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3388,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3390,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3392,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3394,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3396,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3398,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3400,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3402,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3404,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3406,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("exp_K_mult_Q_mult_V", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3408,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3410,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3412,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3414,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3416,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3418,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3420,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3422,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3424,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3426,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3428,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3430,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3432,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3434,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3436,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3438,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("attention", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3440,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3442,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3444,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3446,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3448,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3450,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3452,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3454,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3456,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3458,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3460,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3462,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3464,0,"[0]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3466,0,"[1]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3468,0,"[2]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3470,0,"[3]",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->declBit(c+3472,0,"valid_result",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBus(c+3537,0,"MULT_CYCLES",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBus(c+1194,0,"counter_q",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBus(c+1195,0,"counter_d",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBit(c+1196,0,"valid_result_q",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+193,0,"doProcess_d",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess_q",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable_q",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1199,0,"weight_load_done",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBus(c+1200,0,"weight_load_counter",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::LOGIC, false,-1, 1,0);
    tracep->pushPrefix("data_in_array", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+1201+i*2,0,"",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("K_matrix_in", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+1,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+5,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+7,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+9,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+11,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+13,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+15,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+17,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+19,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+21,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+23,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+25,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+27,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+29,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+31,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("Q_mult_K_result", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+1209+i*2,0,"",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("K_mult_Q_systolicArray", VerilatedTracePrefixType::SCOPE_MODULE);
    VattentionSystolicArray___024root__trace_init_sub__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0(vlSelf, tracep);
    tracep->popPrefix();
    tracep->pushPrefix("exp_Q_mult_K_result", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+1217+i*2,0,"",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->declDouble(c+1225,0,"softmax_norm",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->pushPrefix("Q_mult_K_result_normed", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+194+i*2,0,"",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("data_in_v", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+1227+i*2,0,"",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("weight_V_matrix", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+33,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+35,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+37,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+39,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+41,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+43,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+45,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+47,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+49,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+51,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+53,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+55,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+57,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+59,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+61,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+63,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("exp_Q_mult_K_mult_V_result", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+1235+i*2,0,"",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("V_mult_systolicArray", VerilatedTracePrefixType::SCOPE_MODULE);
    VattentionSystolicArray___024root__trace_init_sub__TOP__attentionSystolicArray__DOT__V_mult_systolicArray__0(vlSelf, tracep);
    tracep->popPrefix();
    tracep->pushPrefix("attention_result", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+1243+i*2,0,"",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("exp_systolic_array", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+3536,0,"K",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBus(c+3535,0,"N",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->pushPrefix("data_in", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+194+i*2,0,"",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("exp_out", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+1217+i*2,0,"",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->declDouble(c+1225,0,"exp_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->pushPrefix("accumulated_sum_wire", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 5; ++i) {
        tracep->declDouble(c+1251+i*2,0,"",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("row[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3473,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+202,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1261,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1263,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1265,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1267,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1269,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+204,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+206,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+208,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[10]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3475,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+210,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1271,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1273,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1275,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1277,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1279,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+212,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+214,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+216,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[11]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3477,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+218,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1281,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1283,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1285,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1287,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1289,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+220,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+222,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+224,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[12]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3479,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+226,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1291,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1293,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1295,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1297,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1299,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+228,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+230,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+232,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[13]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3481,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+234,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1301,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1303,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1305,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1307,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1309,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+236,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+238,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+240,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[14]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3483,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+242,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1311,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1313,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1315,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1317,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1319,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+244,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+246,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+248,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[15]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3485,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+250,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1321,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1323,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1325,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1327,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1329,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+252,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+254,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+256,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[16]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3487,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+258,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1331,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1333,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1335,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1337,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1339,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+260,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+262,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+264,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[17]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3489,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+266,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1341,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1343,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1345,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1347,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1349,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+268,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+270,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+272,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[18]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3491,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+274,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1351,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1353,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1355,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1357,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1359,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+276,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+278,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+280,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[19]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3493,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+282,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1361,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1363,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1365,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1367,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1369,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+284,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+286,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+288,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3495,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+290,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1371,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1373,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1375,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1377,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1379,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+292,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+294,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+296,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[20]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3497,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+298,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1381,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1383,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1385,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1387,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1389,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+300,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+302,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+304,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[21]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3499,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+306,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1391,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1393,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1395,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1397,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1399,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+308,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+310,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+312,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[22]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3501,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+314,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1401,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1403,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1405,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1407,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1409,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+316,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+318,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+320,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[23]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3503,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+322,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1411,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1413,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1415,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1417,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1419,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+324,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+326,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+328,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[24]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3505,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+330,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1421,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1423,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1425,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1427,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1429,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+332,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+334,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+336,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[25]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3507,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+338,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1431,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1433,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1435,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1437,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1439,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+340,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+342,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+344,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[26]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3509,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+346,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1441,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1443,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1445,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1447,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1449,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+348,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+350,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+352,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[27]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3511,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+354,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1451,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1453,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1455,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1457,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1459,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+356,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+358,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+360,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[28]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3513,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+362,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1461,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1463,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1465,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1467,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1469,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+364,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+366,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+368,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[29]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3515,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+370,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1471,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1473,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1475,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1477,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1479,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+372,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+374,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+376,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3517,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+378,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1481,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1483,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1485,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1487,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1489,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+380,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+382,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+384,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[30]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3519,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+386,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1491,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1493,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1495,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1497,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1499,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+388,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+390,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+392,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3521,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+394,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1501,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1503,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1505,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1507,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1509,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+396,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+398,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+400,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[4]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3523,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+402,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1511,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1513,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1515,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1517,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1519,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+404,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+406,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+408,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[5]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3525,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+410,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1521,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1523,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1525,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1527,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1529,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+412,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+414,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+416,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[6]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3527,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+418,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1531,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1533,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1535,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1537,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1539,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+420,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+422,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+424,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[7]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3529,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+426,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1541,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1543,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1545,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1547,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1549,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+428,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+430,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+432,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[8]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3531,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+434,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1551,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1553,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1555,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1557,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1559,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+436,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+438,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+440,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[9]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3533,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+442,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1561,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1563,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1565,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1567,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1569,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+444,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+446,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+448,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("pe_acc", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+1571,0,"exp_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1573,0,"exp_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1575,0,"exp_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1577,0,"exp_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1579,0,"cumulative_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3473,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+450,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1581,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1583,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1585,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1587,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1589,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+452,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+454,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+456,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[10]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3475,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+458,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1591,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1593,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1595,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1597,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1599,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+460,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+462,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+464,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[11]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3477,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+466,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1601,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1603,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1605,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1607,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1609,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+468,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+470,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+472,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[12]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3479,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+474,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1611,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1613,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1615,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1617,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1619,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+476,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+478,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+480,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[13]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3481,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+482,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1621,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1623,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1625,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1627,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1629,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+484,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+486,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+488,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[14]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3483,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+490,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1631,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1633,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1635,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1637,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1639,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+492,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+494,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+496,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[15]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3485,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+498,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1641,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1643,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1645,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1647,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1649,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+500,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+502,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+504,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[16]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3487,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+506,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1651,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1653,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1655,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1657,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1659,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+508,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+510,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+512,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[17]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3489,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+514,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1661,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1663,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1665,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1667,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1669,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+516,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+518,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+520,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[18]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3491,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+522,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1671,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1673,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1675,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1677,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1679,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+524,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+526,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+528,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[19]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3493,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+530,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1681,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1683,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1685,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1687,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1689,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+532,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+534,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+536,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3495,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+538,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1691,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1693,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1695,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1697,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1699,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+540,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+542,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+544,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[20]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3497,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+546,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1701,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1703,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1705,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1707,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1709,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+548,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+550,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+552,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[21]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3499,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+554,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1711,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1713,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1715,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1717,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1719,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+556,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+558,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+560,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[22]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3501,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+562,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1721,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1723,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1725,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1727,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1729,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+564,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+566,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+568,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[23]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3503,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+570,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1731,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1733,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1735,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1737,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1739,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+572,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+574,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+576,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[24]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3505,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+578,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1741,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1743,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1745,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1747,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1749,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+580,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+582,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+584,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[25]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3507,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+586,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1751,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1753,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1755,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1757,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1759,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+588,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+590,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+592,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[26]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3509,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+594,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1761,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1763,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1765,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1767,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1769,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+596,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+598,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+600,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[27]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3511,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+602,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1771,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1773,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1775,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1777,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1779,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+604,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+606,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+608,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[28]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3513,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+610,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1781,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1783,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1785,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1787,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1789,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+612,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+614,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+616,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[29]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3515,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+618,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1791,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1793,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1795,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1797,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1799,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+620,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+622,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+624,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3517,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+626,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1801,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1803,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1805,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1807,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1809,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+628,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+630,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+632,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[30]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3519,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+634,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1811,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1813,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1815,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1817,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1819,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+636,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+638,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+640,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3521,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+642,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1821,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1823,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1825,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1827,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1829,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+644,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+646,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+648,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[4]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3523,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+650,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1831,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1833,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1835,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1837,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1839,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+652,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+654,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+656,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[5]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3525,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+658,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1841,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1843,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1845,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1847,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1849,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+660,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+662,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+664,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[6]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3527,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+666,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1851,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1853,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1855,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1857,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1859,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+668,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+670,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+672,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[7]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3529,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+674,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1861,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1863,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1865,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1867,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1869,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+676,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+678,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+680,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[8]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3531,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+682,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1871,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1873,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1875,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1877,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1879,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+684,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+686,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+688,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[9]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3533,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+690,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1881,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1883,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1885,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1887,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1889,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+692,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+694,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+696,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("pe_acc", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+1891,0,"exp_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1893,0,"exp_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1895,0,"exp_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1897,0,"exp_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1899,0,"cumulative_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3473,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+698,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1901,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1903,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1905,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1907,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1909,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+700,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+702,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+704,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[10]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3475,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+706,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1911,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1913,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1915,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1917,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1919,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+708,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+710,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+712,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[11]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3477,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+714,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1921,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1923,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1925,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1927,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1929,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+716,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+718,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+720,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[12]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3479,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+722,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1931,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1933,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1935,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1937,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1939,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+724,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+726,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+728,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[13]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3481,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+730,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1941,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1943,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1945,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1947,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1949,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+732,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+734,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+736,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[14]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3483,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+738,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1951,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1953,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1955,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1957,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1959,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+740,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+742,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+744,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[15]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3485,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+746,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1961,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1963,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1965,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1967,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1969,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+748,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+750,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+752,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[16]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3487,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+754,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1971,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1973,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1975,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1977,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1979,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+756,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+758,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+760,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[17]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3489,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+762,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1981,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1983,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1985,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1987,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1989,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+764,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+766,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+768,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[18]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3491,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+770,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1991,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1993,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1995,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1997,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1999,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+772,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+774,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+776,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[19]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3493,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+778,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2001,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2003,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2005,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2007,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2009,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+780,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+782,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+784,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3495,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+786,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2011,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2013,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2015,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2017,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2019,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+788,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+790,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+792,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[20]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3497,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+794,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2021,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2023,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2025,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2027,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2029,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+796,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+798,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+800,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[21]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3499,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+802,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2031,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2033,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2035,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2037,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2039,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+804,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+806,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+808,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[22]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3501,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+810,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2041,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2043,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2045,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2047,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2049,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+812,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+814,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+816,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[23]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3503,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+818,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2051,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2053,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2055,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2057,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2059,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+820,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+822,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+824,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[24]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3505,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+826,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2061,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2063,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2065,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2067,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2069,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+828,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+830,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+832,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[25]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3507,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+834,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2071,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2073,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2075,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2077,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2079,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+836,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+838,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+840,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[26]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3509,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+842,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2081,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2083,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2085,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2087,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2089,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+844,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+846,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+848,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[27]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3511,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+850,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2091,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2093,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2095,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2097,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2099,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+852,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+854,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+856,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[28]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3513,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+858,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2101,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2103,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2105,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2107,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2109,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+860,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+862,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+864,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[29]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3515,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+866,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2111,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2113,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2115,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2117,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2119,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+868,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+870,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+872,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3517,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+874,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2121,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2123,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2125,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2127,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2129,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+876,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+878,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+880,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[30]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3519,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+882,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2131,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2133,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2135,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2137,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2139,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+884,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+886,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+888,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3521,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+890,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2141,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2143,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2145,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2147,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2149,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+892,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+894,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+896,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[4]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3523,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+898,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2151,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2153,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2155,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2157,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2159,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+900,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+902,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+904,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[5]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3525,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+906,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2161,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2163,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2165,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2167,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2169,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+908,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+910,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+912,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[6]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3527,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+914,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2171,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2173,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2175,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2177,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2179,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+916,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+918,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+920,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[7]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3529,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+922,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2181,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2183,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2185,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2187,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2189,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+924,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+926,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+928,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[8]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3531,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+930,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2191,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2193,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2195,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2197,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2199,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+932,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+934,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+936,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[9]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3533,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+938,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2201,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2203,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2205,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2207,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2209,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+940,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+942,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+944,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("pe_acc", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+2211,0,"exp_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2213,0,"exp_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2215,0,"exp_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2217,0,"exp_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2219,0,"cumulative_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3473,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+946,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2221,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2223,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2225,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2227,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2229,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+948,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+950,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+952,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[10]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3475,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+954,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2231,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2233,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2235,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2237,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2239,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+956,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+958,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+960,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[11]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3477,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+962,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2241,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2243,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2245,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2247,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2249,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+964,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+966,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+968,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[12]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3479,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+970,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2251,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2253,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2255,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2257,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2259,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+972,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+974,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+976,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[13]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3481,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+978,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2261,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2263,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2265,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2267,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2269,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+980,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+982,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+984,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[14]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3483,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+986,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2271,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2273,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2275,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2277,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2279,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+988,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+990,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+992,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[15]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3485,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+994,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2281,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2283,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2285,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2287,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2289,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+996,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+998,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1000,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[16]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3487,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1002,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2291,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2293,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2295,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2297,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2299,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1004,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1006,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1008,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[17]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3489,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1010,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2301,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2303,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2305,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2307,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2309,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1012,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1014,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1016,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[18]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3491,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1018,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2311,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2313,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2315,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2317,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2319,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1020,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1022,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1024,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[19]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3493,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1026,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2321,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2323,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2325,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2327,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2329,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1028,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1030,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1032,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3495,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1034,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2331,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2333,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2335,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2337,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2339,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1036,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1038,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1040,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[20]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3497,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1042,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2341,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2343,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2345,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2347,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2349,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1044,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1046,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1048,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[21]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3499,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1050,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2351,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2353,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2355,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2357,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2359,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1052,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1054,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1056,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[22]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3501,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1058,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2361,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2363,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2365,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2367,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2369,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1060,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1062,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1064,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[23]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3503,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1066,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2371,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2373,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2375,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2377,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2379,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1068,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1070,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1072,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[24]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3505,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1074,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2381,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2383,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2385,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2387,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2389,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1076,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1078,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1080,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[25]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3507,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1082,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2391,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2393,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2395,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2397,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2399,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1084,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1086,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1088,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[26]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3509,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1090,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2401,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2403,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2405,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2407,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2409,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1092,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1094,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1096,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[27]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3511,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1098,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2411,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2413,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2415,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2417,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2419,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1100,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1102,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1104,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[28]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3513,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1106,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2421,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2423,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2425,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2427,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2429,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1108,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1110,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1112,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[29]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3515,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1114,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2431,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2433,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2435,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2437,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2439,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1116,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1118,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1120,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3517,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1122,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2441,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2443,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2445,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2447,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2449,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1124,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1126,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1128,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[30]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3519,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1130,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2451,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2453,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2455,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2457,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2459,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1132,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1134,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1136,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3521,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1138,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2461,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2463,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2465,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2467,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2469,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1140,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1142,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1144,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[4]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3523,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1146,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2471,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2473,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2475,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2477,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2479,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1148,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1150,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1152,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[5]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3525,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1154,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2481,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2483,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2485,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2487,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2489,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1156,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1158,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1160,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[6]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3527,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1162,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2491,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2493,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2495,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2497,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2499,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1164,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1166,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1168,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[7]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3529,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1170,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2501,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2503,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2505,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2507,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2509,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1172,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1174,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1176,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[8]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3531,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1178,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2511,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2513,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2515,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2517,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2519,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1180,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1182,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1184,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[9]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+3533,0,"k_factorial",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1186,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2521,0,"x_pow_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2523,0,"taylor_approx_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2525,0,"x_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2527,0,"x_pow_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2529,0,"taylor_approx_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1188,0,"exp_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1190,0,"taylor_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+1192,0,"taylor_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("pe_acc", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+2531,0,"exp_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2533,0,"exp_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2535,0,"exp_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2537,0,"exp_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2539,0,"cumulative_sum",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("softmax_norm_systolic_array", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+3535,0,"N",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+1225,0,"norm_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->pushPrefix("data_in", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+1235+i*2,0,"",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("data_normed_out", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+1243+i*2,0,"",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("norm_wire", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 5; ++i) {
        tracep->declDouble(c+2541+i*2,0,"",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("dummy_pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+1225,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2551,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("row[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+2553,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2555,0,"norm_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2557,0,"x_normed_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2559,0,"norm_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2561,0,"normed_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+2563,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2565,0,"norm_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2567,0,"x_normed_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2569,0,"norm_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2571,0,"normed_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+2573,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2575,0,"norm_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2577,0,"x_normed_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2579,0,"norm_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2581,0,"normed_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+2583,0,"x_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2585,0,"norm_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2587,0,"x_normed_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2589,0,"norm_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2591,0,"normed_term",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("unnamedblk1", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+2593,0,"j",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->popPrefix();
    tracep->pushPrefix("unnamedblk11", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+2594,0,"i",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->pushPrefix("unnamedblk12", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+2595,0,"j",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("unnamedblk13", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+2596,0,"i",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->pushPrefix("unnamedblk14", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+2597,0,"j",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("unnamedblk15", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+2598,0,"i",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->pushPrefix("unnamedblk16", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+2599,0,"j",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("unnamedblk2", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+2600,0,"j",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->popPrefix();
    tracep->pushPrefix("unnamedblk3", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+3538,0,"i",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->pushPrefix("unnamedblk4", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+3538,0,"j",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("unnamedblk5", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+3538,0,"i",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->popPrefix();
    tracep->pushPrefix("unnamedblk6", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+3538,0,"i",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->popPrefix();
    tracep->pushPrefix("unnamedblk7", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+3538,0,"i",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->pushPrefix("unnamedblk8", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+3538,0,"j",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("unnamedblk9", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+2601,0,"i",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->pushPrefix("unnamedblk10", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBus(c+2602,0,"j",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
}

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_init_sub__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0(VattentionSystolicArray___024root* vlSelf, VerilatedVcd* tracep) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root__trace_init_sub__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    const int c = vlSymsp->__Vm_baseCode;
    // Body
    tracep->declBus(c+3535,0,"N",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBus(c+3535,0,"M",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->pushPrefix("weight_in", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+65,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+67,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+69,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+71,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+73,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+75,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+77,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+79,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+81,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+83,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+85,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+87,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+89,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+91,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+93,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+95,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("data_in", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+2603+i*2,0,"",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("result_out", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+2611+i*2,0,"",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("data_wire", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2619,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2621,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2623,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2625,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2627,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2629,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2631,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2633,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2635,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2637,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2639,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2641,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2643,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2645,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2647,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2649,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[4]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2651,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2653,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2655,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2657,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("partial_sum_wire", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2659,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2661,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2663,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2665,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2667,0,"[4]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2669,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2671,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2673,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2675,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2677,0,"[4]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2679,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2681,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2683,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2685,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2687,0,"[4]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2689,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2691,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2693,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2695,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2697,0,"[4]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+97,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2699,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2701,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2703,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2705,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2707,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2709,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2711,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+99,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2713,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2715,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2717,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2719,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2721,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2723,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2725,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+101,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2727,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2729,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2731,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2733,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2735,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2737,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2739,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+103,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2741,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2743,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2745,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2747,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2749,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2751,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2753,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+105,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2755,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2757,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2759,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2761,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2763,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2765,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2767,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+107,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2769,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2771,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2773,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2775,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2777,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2779,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2781,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+109,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2783,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2785,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2787,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2789,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2791,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2793,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2795,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+111,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2797,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2799,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2801,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2803,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2805,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2807,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2809,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+113,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2811,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2813,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2815,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2817,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2819,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2821,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2823,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+115,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2825,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2827,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2829,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2831,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2833,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2835,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2837,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+117,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2839,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2841,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2843,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2845,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2847,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2849,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2851,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+119,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2853,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2855,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2857,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2859,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2861,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2863,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2865,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+121,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2867,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2869,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2871,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2873,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2875,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2877,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2879,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+123,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2881,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2883,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2885,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2887,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2889,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2891,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2893,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+125,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2895,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2897,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2899,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2901,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2903,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2905,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2907,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+127,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2909,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2911,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2913,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2915,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2917,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2919,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2921,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
}

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_init_sub__TOP__attentionSystolicArray__DOT__V_mult_systolicArray__0(VattentionSystolicArray___024root* vlSelf, VerilatedVcd* tracep) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root__trace_init_sub__TOP__attentionSystolicArray__DOT__V_mult_systolicArray__0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    const int c = vlSymsp->__Vm_baseCode;
    // Body
    tracep->declBus(c+3535,0,"N",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBus(c+3535,0,"M",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::PARAMETER, VerilatedTraceSigType::INT, false,-1, 31,0);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->pushPrefix("weight_in", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+129,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+131,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+133,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+135,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+137,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+139,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+141,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+143,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+145,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+147,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+149,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+151,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+153,0,"[0]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+155,0,"[1]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+157,0,"[2]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+159,0,"[3]",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("data_in", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+2923+i*2,0,"",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("result_out", VerilatedTracePrefixType::ARRAY_UNPACKED);
    for (int i = 0; i < 4; ++i) {
        tracep->declDouble(c+2931+i*2,0,"",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, true,(i+0));
    }
    tracep->popPrefix();
    tracep->pushPrefix("data_wire", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2939,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2941,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2943,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2945,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2947,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2949,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2951,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2953,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2955,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2957,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2959,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2961,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2963,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2965,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2967,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2969,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[4]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2971,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2973,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2975,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2977,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("partial_sum_wire", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->pushPrefix("[0]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2979,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2981,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2983,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2985,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2987,0,"[4]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[1]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2989,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2991,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2993,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2995,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+2997,0,"[4]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[2]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+2999,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3001,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3003,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3005,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3007,0,"[4]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->pushPrefix("[3]", VerilatedTracePrefixType::ARRAY_UNPACKED);
    tracep->declDouble(c+3009,0,"[0]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3011,0,"[1]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3013,0,"[2]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3015,0,"[3]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3017,0,"[4]",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+161,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3019,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3021,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3023,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3025,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3027,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3029,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3031,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+163,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3033,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3035,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3037,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3039,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3041,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3043,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3045,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+165,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3047,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3049,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3051,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3053,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3055,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3057,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3059,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+167,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3061,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3063,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3065,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3067,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3069,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3071,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3073,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+169,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3075,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3077,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3079,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3081,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3083,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3085,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3087,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+171,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3089,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3091,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3093,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3095,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3097,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3099,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3101,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+173,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3103,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3105,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3107,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3109,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3111,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3113,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3115,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+175,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3117,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3119,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3121,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3123,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3125,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3127,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3129,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+177,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3131,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3133,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3135,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3137,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3139,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3141,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3143,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+179,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3145,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3147,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3149,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3151,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3153,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3155,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3157,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+181,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3159,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3161,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3163,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3165,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3167,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3169,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3171,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+183,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3173,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3175,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3177,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3179,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3181,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3183,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3185,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("row[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("col[0]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+185,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3187,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3189,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3191,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3193,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3195,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3197,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3199,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[1]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+187,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3201,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3203,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3205,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3207,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3209,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3211,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3213,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[2]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+189,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3215,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3217,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3219,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3221,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3223,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3225,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3227,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->pushPrefix("col[3]", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->pushPrefix("pe", VerilatedTracePrefixType::SCOPE_MODULE);
    tracep->declBit(c+3243,0,"clk",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+3244,0,"reset",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1197,0,"doProcess",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declBit(c+1198,0,"weight_load_enable",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::LOGIC, false,-1);
    tracep->declDouble(c+191,0,"weight_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3229,0,"data_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3231,0,"partial_sum_in",-1, VerilatedTraceSigDirection::INPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3233,0,"data_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3235,0,"partial_sum_out",-1, VerilatedTraceSigDirection::OUTPUT, VerilatedTraceSigKind::WIRE, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3237,0,"weight",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3239,0,"mul_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->declDouble(c+3241,0,"add_result",-1, VerilatedTraceSigDirection::NONE, VerilatedTraceSigKind::VAR, VerilatedTraceSigType::DOUBLE, false,-1);
    tracep->popPrefix();
    tracep->popPrefix();
    tracep->popPrefix();
}

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_init_top(VattentionSystolicArray___024root* vlSelf, VerilatedVcd* tracep) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root__trace_init_top\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    VattentionSystolicArray___024root__trace_init_sub__TOP__0(vlSelf, tracep);
}

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_const_0(void* voidSelf, VerilatedVcd::Buffer* bufp);
VL_ATTR_COLD void VattentionSystolicArray___024root__trace_full_0(void* voidSelf, VerilatedVcd::Buffer* bufp);
void VattentionSystolicArray___024root__trace_chg_0(void* voidSelf, VerilatedVcd::Buffer* bufp);
void VattentionSystolicArray___024root__trace_cleanup(void* voidSelf, VerilatedVcd* /*unused*/);

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_register(VattentionSystolicArray___024root* vlSelf, VerilatedVcd* tracep) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root__trace_register\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    tracep->addConstCb(&VattentionSystolicArray___024root__trace_const_0, 0U, vlSelf);
    tracep->addFullCb(&VattentionSystolicArray___024root__trace_full_0, 0U, vlSelf);
    tracep->addChgCb(&VattentionSystolicArray___024root__trace_chg_0, 0U, vlSelf);
    tracep->addCleanupCb(&VattentionSystolicArray___024root__trace_cleanup, vlSelf);
}

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_const_0_sub_0(VattentionSystolicArray___024root* vlSelf, VerilatedVcd::Buffer* bufp);

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_const_0(void* voidSelf, VerilatedVcd::Buffer* bufp) {
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root__trace_const_0\n"); );
    // Init
    VattentionSystolicArray___024root* const __restrict vlSelf VL_ATTR_UNUSED = static_cast<VattentionSystolicArray___024root*>(voidSelf);
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    // Body
    VattentionSystolicArray___024root__trace_const_0_sub_0((&vlSymsp->TOP), bufp);
}

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_const_0_sub_0(VattentionSystolicArray___024root* vlSelf, VerilatedVcd::Buffer* bufp) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root__trace_const_0_sub_0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    uint32_t* const oldp VL_ATTR_UNUSED = bufp->oldp(vlSymsp->__Vm_baseCode);
    // Body
    bufp->fullIData(oldp+3535,(4U),32);
    bufp->fullIData(oldp+3536,(0x20U),32);
    bufp->fullIData(oldp+3537,(0x31U),32);
    bufp->fullIData(oldp+3538,(4U),32);
}

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_full_0_sub_0(VattentionSystolicArray___024root* vlSelf, VerilatedVcd::Buffer* bufp);

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_full_0(void* voidSelf, VerilatedVcd::Buffer* bufp) {
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root__trace_full_0\n"); );
    // Init
    VattentionSystolicArray___024root* const __restrict vlSelf VL_ATTR_UNUSED = static_cast<VattentionSystolicArray___024root*>(voidSelf);
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    // Body
    VattentionSystolicArray___024root__trace_full_0_sub_0((&vlSymsp->TOP), bufp);
}

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_full_0_sub_0(VattentionSystolicArray___024root* vlSelf, VerilatedVcd::Buffer* bufp) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root__trace_full_0_sub_0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    uint32_t* const oldp VL_ATTR_UNUSED = bufp->oldp(vlSymsp->__Vm_baseCode);
    // Body
    bufp->fullDouble(oldp+1,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                             [0U][0U]));
    bufp->fullDouble(oldp+3,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                             [0U][1U]));
    bufp->fullDouble(oldp+5,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                             [0U][2U]));
    bufp->fullDouble(oldp+7,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                             [0U][3U]));
    bufp->fullDouble(oldp+9,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                             [1U][0U]));
    bufp->fullDouble(oldp+11,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [1U][1U]));
    bufp->fullDouble(oldp+13,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [1U][2U]));
    bufp->fullDouble(oldp+15,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [1U][3U]));
    bufp->fullDouble(oldp+17,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [2U][0U]));
    bufp->fullDouble(oldp+19,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [2U][1U]));
    bufp->fullDouble(oldp+21,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [2U][2U]));
    bufp->fullDouble(oldp+23,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [2U][3U]));
    bufp->fullDouble(oldp+25,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [3U][0U]));
    bufp->fullDouble(oldp+27,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [3U][1U]));
    bufp->fullDouble(oldp+29,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [3U][2U]));
    bufp->fullDouble(oldp+31,(vlSelfRef.attentionSystolicArray__DOT__K_matrix_in
                              [3U][3U]));
    bufp->fullDouble(oldp+33,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [0U][0U]));
    bufp->fullDouble(oldp+35,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [0U][1U]));
    bufp->fullDouble(oldp+37,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [0U][2U]));
    bufp->fullDouble(oldp+39,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [0U][3U]));
    bufp->fullDouble(oldp+41,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [1U][0U]));
    bufp->fullDouble(oldp+43,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [1U][1U]));
    bufp->fullDouble(oldp+45,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [1U][2U]));
    bufp->fullDouble(oldp+47,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [1U][3U]));
    bufp->fullDouble(oldp+49,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [2U][0U]));
    bufp->fullDouble(oldp+51,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [2U][1U]));
    bufp->fullDouble(oldp+53,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [2U][2U]));
    bufp->fullDouble(oldp+55,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [2U][3U]));
    bufp->fullDouble(oldp+57,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [3U][0U]));
    bufp->fullDouble(oldp+59,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [3U][1U]));
    bufp->fullDouble(oldp+61,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [3U][2U]));
    bufp->fullDouble(oldp+63,(vlSelfRef.attentionSystolicArray__DOT__weight_V_matrix
                              [3U][3U]));
    bufp->fullDouble(oldp+65,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [0U][0U]));
    bufp->fullDouble(oldp+67,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [0U][1U]));
    bufp->fullDouble(oldp+69,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [0U][2U]));
    bufp->fullDouble(oldp+71,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [0U][3U]));
    bufp->fullDouble(oldp+73,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [1U][0U]));
    bufp->fullDouble(oldp+75,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [1U][1U]));
    bufp->fullDouble(oldp+77,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [1U][2U]));
    bufp->fullDouble(oldp+79,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [1U][3U]));
    bufp->fullDouble(oldp+81,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [2U][0U]));
    bufp->fullDouble(oldp+83,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [2U][1U]));
    bufp->fullDouble(oldp+85,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [2U][2U]));
    bufp->fullDouble(oldp+87,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [2U][3U]));
    bufp->fullDouble(oldp+89,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [3U][0U]));
    bufp->fullDouble(oldp+91,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [3U][1U]));
    bufp->fullDouble(oldp+93,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [3U][2U]));
    bufp->fullDouble(oldp+95,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [3U][3U]));
    bufp->fullDouble(oldp+97,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [0U][0U]));
    bufp->fullDouble(oldp+99,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                              [0U][1U]));
    bufp->fullDouble(oldp+101,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [0U][2U]));
    bufp->fullDouble(oldp+103,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [0U][3U]));
    bufp->fullDouble(oldp+105,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [1U][0U]));
    bufp->fullDouble(oldp+107,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [1U][1U]));
    bufp->fullDouble(oldp+109,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [1U][2U]));
    bufp->fullDouble(oldp+111,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [1U][3U]));
    bufp->fullDouble(oldp+113,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [2U][0U]));
    bufp->fullDouble(oldp+115,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [2U][1U]));
    bufp->fullDouble(oldp+117,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [2U][2U]));
    bufp->fullDouble(oldp+119,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [2U][3U]));
    bufp->fullDouble(oldp+121,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [3U][0U]));
    bufp->fullDouble(oldp+123,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [3U][1U]));
    bufp->fullDouble(oldp+125,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [3U][2U]));
    bufp->fullDouble(oldp+127,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.weight_in
                               [3U][3U]));
    bufp->fullDouble(oldp+129,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [0U][0U]));
    bufp->fullDouble(oldp+131,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [0U][1U]));
    bufp->fullDouble(oldp+133,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [0U][2U]));
    bufp->fullDouble(oldp+135,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [0U][3U]));
    bufp->fullDouble(oldp+137,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [1U][0U]));
    bufp->fullDouble(oldp+139,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [1U][1U]));
    bufp->fullDouble(oldp+141,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [1U][2U]));
    bufp->fullDouble(oldp+143,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [1U][3U]));
    bufp->fullDouble(oldp+145,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [2U][0U]));
    bufp->fullDouble(oldp+147,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [2U][1U]));
    bufp->fullDouble(oldp+149,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [2U][2U]));
    bufp->fullDouble(oldp+151,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [2U][3U]));
    bufp->fullDouble(oldp+153,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [3U][0U]));
    bufp->fullDouble(oldp+155,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [3U][1U]));
    bufp->fullDouble(oldp+157,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [3U][2U]));
    bufp->fullDouble(oldp+159,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [3U][3U]));
    bufp->fullDouble(oldp+161,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [0U][0U]));
    bufp->fullDouble(oldp+163,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [0U][1U]));
    bufp->fullDouble(oldp+165,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [0U][2U]));
    bufp->fullDouble(oldp+167,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [0U][3U]));
    bufp->fullDouble(oldp+169,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [1U][0U]));
    bufp->fullDouble(oldp+171,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [1U][1U]));
    bufp->fullDouble(oldp+173,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [1U][2U]));
    bufp->fullDouble(oldp+175,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [1U][3U]));
    bufp->fullDouble(oldp+177,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [2U][0U]));
    bufp->fullDouble(oldp+179,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [2U][1U]));
    bufp->fullDouble(oldp+181,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [2U][2U]));
    bufp->fullDouble(oldp+183,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [2U][3U]));
    bufp->fullDouble(oldp+185,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [3U][0U]));
    bufp->fullDouble(oldp+187,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [3U][1U]));
    bufp->fullDouble(oldp+189,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [3U][2U]));
    bufp->fullDouble(oldp+191,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.weight_in
                               [3U][3U]));
    bufp->fullBit(oldp+193,(vlSelfRef.attentionSystolicArray__DOT__doProcess_d));
    bufp->fullDouble(oldp+194,(vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed[0]));
    bufp->fullDouble(oldp+196,(vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed[1]));
    bufp->fullDouble(oldp+198,(vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed[2]));
    bufp->fullDouble(oldp+200,(vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result_normed[3]));
    bufp->fullDouble(oldp+202,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0U]));
    bufp->fullDouble(oldp+204,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+206,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+208,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+210,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0xaU]));
    bufp->fullDouble(oldp+212,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+214,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+216,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+218,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0xbU]));
    bufp->fullDouble(oldp+220,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+222,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+224,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+226,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0xcU]));
    bufp->fullDouble(oldp+228,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+230,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+232,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+234,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0xdU]));
    bufp->fullDouble(oldp+236,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+238,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+240,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+242,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0xeU]));
    bufp->fullDouble(oldp+244,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+246,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+248,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+250,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0xfU]));
    bufp->fullDouble(oldp+252,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+254,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+256,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+258,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x10U]));
    bufp->fullDouble(oldp+260,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+262,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+264,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+266,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x11U]));
    bufp->fullDouble(oldp+268,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+270,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+272,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+274,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x12U]));
    bufp->fullDouble(oldp+276,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+278,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+280,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+282,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x13U]));
    bufp->fullDouble(oldp+284,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+286,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+288,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+290,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][1U]));
    bufp->fullDouble(oldp+292,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+294,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+296,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+298,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x14U]));
    bufp->fullDouble(oldp+300,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+302,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+304,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+306,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x15U]));
    bufp->fullDouble(oldp+308,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+310,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+312,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+314,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x16U]));
    bufp->fullDouble(oldp+316,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+318,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+320,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+322,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x17U]));
    bufp->fullDouble(oldp+324,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+326,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+328,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+330,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x18U]));
    bufp->fullDouble(oldp+332,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+334,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+336,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+338,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x19U]));
    bufp->fullDouble(oldp+340,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+342,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+344,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+346,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x1aU]));
    bufp->fullDouble(oldp+348,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+350,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+352,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+354,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x1bU]));
    bufp->fullDouble(oldp+356,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+358,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+360,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+362,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x1cU]));
    bufp->fullDouble(oldp+364,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+366,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+368,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+370,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x1dU]));
    bufp->fullDouble(oldp+372,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+374,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+376,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+378,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][2U]));
    bufp->fullDouble(oldp+380,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+382,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+384,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+386,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][0x1eU]));
    bufp->fullDouble(oldp+388,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+390,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+392,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+394,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][3U]));
    bufp->fullDouble(oldp+396,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+398,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+400,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+402,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][4U]));
    bufp->fullDouble(oldp+404,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+406,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+408,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+410,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][5U]));
    bufp->fullDouble(oldp+412,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+414,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+416,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+418,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][6U]));
    bufp->fullDouble(oldp+420,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+422,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+424,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+426,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][7U]));
    bufp->fullDouble(oldp+428,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+430,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+432,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+434,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][8U]));
    bufp->fullDouble(oldp+436,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+438,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+440,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+442,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [0U][9U]));
    bufp->fullDouble(oldp+444,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+446,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+448,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+450,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0U]));
    bufp->fullDouble(oldp+452,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+454,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+456,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+458,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0xaU]));
    bufp->fullDouble(oldp+460,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+462,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+464,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+466,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0xbU]));
    bufp->fullDouble(oldp+468,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+470,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+472,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+474,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0xcU]));
    bufp->fullDouble(oldp+476,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+478,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+480,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+482,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0xdU]));
    bufp->fullDouble(oldp+484,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+486,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+488,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+490,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0xeU]));
    bufp->fullDouble(oldp+492,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+494,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+496,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+498,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0xfU]));
    bufp->fullDouble(oldp+500,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+502,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+504,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+506,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x10U]));
    bufp->fullDouble(oldp+508,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+510,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+512,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+514,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x11U]));
    bufp->fullDouble(oldp+516,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+518,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+520,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+522,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x12U]));
    bufp->fullDouble(oldp+524,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+526,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+528,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+530,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x13U]));
    bufp->fullDouble(oldp+532,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+534,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+536,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+538,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][1U]));
    bufp->fullDouble(oldp+540,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+542,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+544,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+546,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x14U]));
    bufp->fullDouble(oldp+548,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+550,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+552,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+554,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x15U]));
    bufp->fullDouble(oldp+556,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+558,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+560,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+562,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x16U]));
    bufp->fullDouble(oldp+564,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+566,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+568,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+570,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x17U]));
    bufp->fullDouble(oldp+572,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+574,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+576,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+578,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x18U]));
    bufp->fullDouble(oldp+580,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+582,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+584,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+586,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x19U]));
    bufp->fullDouble(oldp+588,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+590,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+592,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+594,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x1aU]));
    bufp->fullDouble(oldp+596,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+598,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+600,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+602,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x1bU]));
    bufp->fullDouble(oldp+604,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+606,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+608,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+610,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x1cU]));
    bufp->fullDouble(oldp+612,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+614,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+616,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+618,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x1dU]));
    bufp->fullDouble(oldp+620,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+622,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+624,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+626,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][2U]));
    bufp->fullDouble(oldp+628,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+630,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+632,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+634,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][0x1eU]));
    bufp->fullDouble(oldp+636,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+638,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+640,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+642,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][3U]));
    bufp->fullDouble(oldp+644,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+646,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+648,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+650,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][4U]));
    bufp->fullDouble(oldp+652,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+654,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+656,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+658,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][5U]));
    bufp->fullDouble(oldp+660,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+662,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+664,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+666,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][6U]));
    bufp->fullDouble(oldp+668,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+670,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+672,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+674,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][7U]));
    bufp->fullDouble(oldp+676,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+678,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+680,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+682,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][8U]));
    bufp->fullDouble(oldp+684,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+686,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+688,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+690,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [1U][9U]));
    bufp->fullDouble(oldp+692,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+694,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+696,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+698,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0U]));
    bufp->fullDouble(oldp+700,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+702,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+704,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+706,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0xaU]));
    bufp->fullDouble(oldp+708,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+710,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+712,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+714,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0xbU]));
    bufp->fullDouble(oldp+716,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+718,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+720,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+722,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0xcU]));
    bufp->fullDouble(oldp+724,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+726,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+728,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+730,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0xdU]));
    bufp->fullDouble(oldp+732,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+734,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+736,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+738,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0xeU]));
    bufp->fullDouble(oldp+740,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+742,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+744,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+746,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0xfU]));
    bufp->fullDouble(oldp+748,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+750,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+752,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+754,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x10U]));
    bufp->fullDouble(oldp+756,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+758,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+760,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+762,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x11U]));
    bufp->fullDouble(oldp+764,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+766,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+768,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+770,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x12U]));
    bufp->fullDouble(oldp+772,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+774,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+776,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+778,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x13U]));
    bufp->fullDouble(oldp+780,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+782,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+784,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+786,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][1U]));
    bufp->fullDouble(oldp+788,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+790,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+792,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+794,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x14U]));
    bufp->fullDouble(oldp+796,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+798,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+800,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+802,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x15U]));
    bufp->fullDouble(oldp+804,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+806,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+808,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+810,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x16U]));
    bufp->fullDouble(oldp+812,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+814,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+816,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+818,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x17U]));
    bufp->fullDouble(oldp+820,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+822,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+824,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+826,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x18U]));
    bufp->fullDouble(oldp+828,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+830,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+832,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+834,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x19U]));
    bufp->fullDouble(oldp+836,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+838,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+840,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+842,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x1aU]));
    bufp->fullDouble(oldp+844,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+846,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+848,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+850,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x1bU]));
    bufp->fullDouble(oldp+852,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+854,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+856,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+858,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x1cU]));
    bufp->fullDouble(oldp+860,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+862,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+864,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+866,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x1dU]));
    bufp->fullDouble(oldp+868,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+870,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+872,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+874,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][2U]));
    bufp->fullDouble(oldp+876,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+878,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+880,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+882,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][0x1eU]));
    bufp->fullDouble(oldp+884,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+886,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+888,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+890,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][3U]));
    bufp->fullDouble(oldp+892,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+894,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+896,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+898,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][4U]));
    bufp->fullDouble(oldp+900,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+902,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+904,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+906,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][5U]));
    bufp->fullDouble(oldp+908,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+910,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+912,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+914,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][6U]));
    bufp->fullDouble(oldp+916,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+918,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+920,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+922,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][7U]));
    bufp->fullDouble(oldp+924,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+926,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+928,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+930,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][8U]));
    bufp->fullDouble(oldp+932,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+934,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+936,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+938,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [2U][9U]));
    bufp->fullDouble(oldp+940,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+942,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+944,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+946,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [3U][0U]));
    bufp->fullDouble(oldp+948,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+950,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+952,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+954,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [3U][0xaU]));
    bufp->fullDouble(oldp+956,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+958,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+960,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+962,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [3U][0xbU]));
    bufp->fullDouble(oldp+964,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+966,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+968,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+970,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [3U][0xcU]));
    bufp->fullDouble(oldp+972,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+974,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+976,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+978,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [3U][0xdU]));
    bufp->fullDouble(oldp+980,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+982,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+984,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+986,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [3U][0xeU]));
    bufp->fullDouble(oldp+988,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+990,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+992,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+994,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                               [3U][0xfU]));
    bufp->fullDouble(oldp+996,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+998,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1000,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1002,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x10U]));
    bufp->fullDouble(oldp+1004,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1006,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1008,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1010,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x11U]));
    bufp->fullDouble(oldp+1012,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1014,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1016,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1018,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x12U]));
    bufp->fullDouble(oldp+1020,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1022,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1024,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1026,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x13U]));
    bufp->fullDouble(oldp+1028,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1030,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1032,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1034,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+1036,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1038,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1040,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1042,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x14U]));
    bufp->fullDouble(oldp+1044,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1046,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1048,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1050,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x15U]));
    bufp->fullDouble(oldp+1052,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1054,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1056,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1058,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x16U]));
    bufp->fullDouble(oldp+1060,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1062,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1064,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1066,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x17U]));
    bufp->fullDouble(oldp+1068,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1070,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1072,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1074,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x18U]));
    bufp->fullDouble(oldp+1076,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1078,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1080,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1082,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x19U]));
    bufp->fullDouble(oldp+1084,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1086,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1088,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1090,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x1aU]));
    bufp->fullDouble(oldp+1092,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1094,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1096,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1098,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x1bU]));
    bufp->fullDouble(oldp+1100,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1102,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1104,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1106,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x1cU]));
    bufp->fullDouble(oldp+1108,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1110,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1112,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1114,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x1dU]));
    bufp->fullDouble(oldp+1116,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1118,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1120,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1122,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+1124,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1126,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1128,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1130,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][0x1eU]));
    bufp->fullDouble(oldp+1132,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1134,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1136,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1138,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+1140,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1142,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1144,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1146,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][4U]));
    bufp->fullDouble(oldp+1148,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1150,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1152,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1154,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][5U]));
    bufp->fullDouble(oldp+1156,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1158,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1160,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1162,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][6U]));
    bufp->fullDouble(oldp+1164,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1166,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1168,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1170,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][7U]));
    bufp->fullDouble(oldp+1172,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1174,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1176,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1178,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][8U]));
    bufp->fullDouble(oldp+1180,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1182,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1184,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullDouble(oldp+1186,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                                [3U][9U]));
    bufp->fullDouble(oldp+1188,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term));
    bufp->fullDouble(oldp+1190,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_term));
    bufp->fullDouble(oldp+1192,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum));
    bufp->fullIData(oldp+1194,(vlSelfRef.attentionSystolicArray__DOT__counter_q),32);
    bufp->fullIData(oldp+1195,(vlSelfRef.attentionSystolicArray__DOT__counter_d),32);
    bufp->fullBit(oldp+1196,(vlSelfRef.attentionSystolicArray__DOT__valid_result_q));
    bufp->fullBit(oldp+1197,(vlSelfRef.attentionSystolicArray__DOT__doProcess_q));
    bufp->fullBit(oldp+1198,(vlSelfRef.attentionSystolicArray__DOT__weight_load_enable_q));
    bufp->fullBit(oldp+1199,((1U & (~ (IData)(vlSelfRef.attentionSystolicArray__DOT__weight_load_enable_q)))));
    bufp->fullCData(oldp+1200,(vlSelfRef.attentionSystolicArray__DOT__weight_load_counter),2);
    bufp->fullDouble(oldp+1201,(vlSelfRef.attentionSystolicArray__DOT__data_in_array[0]));
    bufp->fullDouble(oldp+1203,(vlSelfRef.attentionSystolicArray__DOT__data_in_array[1]));
    bufp->fullDouble(oldp+1205,(vlSelfRef.attentionSystolicArray__DOT__data_in_array[2]));
    bufp->fullDouble(oldp+1207,(vlSelfRef.attentionSystolicArray__DOT__data_in_array[3]));
    bufp->fullDouble(oldp+1209,(vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result[0]));
    bufp->fullDouble(oldp+1211,(vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result[1]));
    bufp->fullDouble(oldp+1213,(vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result[2]));
    bufp->fullDouble(oldp+1215,(vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result[3]));
    bufp->fullDouble(oldp+1217,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result[0]));
    bufp->fullDouble(oldp+1219,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result[1]));
    bufp->fullDouble(oldp+1221,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result[2]));
    bufp->fullDouble(oldp+1223,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result[3]));
    bufp->fullDouble(oldp+1225,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                                [0U]));
    bufp->fullDouble(oldp+1227,(vlSelfRef.attentionSystolicArray__DOT__data_in_v[0]));
    bufp->fullDouble(oldp+1229,(vlSelfRef.attentionSystolicArray__DOT__data_in_v[1]));
    bufp->fullDouble(oldp+1231,(vlSelfRef.attentionSystolicArray__DOT__data_in_v[2]));
    bufp->fullDouble(oldp+1233,(vlSelfRef.attentionSystolicArray__DOT__data_in_v[3]));
    bufp->fullDouble(oldp+1235,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result[0]));
    bufp->fullDouble(oldp+1237,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result[1]));
    bufp->fullDouble(oldp+1239,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result[2]));
    bufp->fullDouble(oldp+1241,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result[3]));
    bufp->fullDouble(oldp+1243,(vlSelfRef.attentionSystolicArray__DOT__attention_result[0]));
    bufp->fullDouble(oldp+1245,(vlSelfRef.attentionSystolicArray__DOT__attention_result[1]));
    bufp->fullDouble(oldp+1247,(vlSelfRef.attentionSystolicArray__DOT__attention_result[2]));
    bufp->fullDouble(oldp+1249,(vlSelfRef.attentionSystolicArray__DOT__attention_result[3]));
    bufp->fullDouble(oldp+1251,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[0]));
    bufp->fullDouble(oldp+1253,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[1]));
    bufp->fullDouble(oldp+1255,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[2]));
    bufp->fullDouble(oldp+1257,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[3]));
    bufp->fullDouble(oldp+1259,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[4]));
    bufp->fullDouble(oldp+1261,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0U]));
    bufp->fullDouble(oldp+1263,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0U]));
    bufp->fullDouble(oldp+1265,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1267,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1269,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1271,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0xaU]));
    bufp->fullDouble(oldp+1273,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0xaU]));
    bufp->fullDouble(oldp+1275,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1277,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1279,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1281,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0xbU]));
    bufp->fullDouble(oldp+1283,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0xbU]));
    bufp->fullDouble(oldp+1285,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1287,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1289,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1291,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0xcU]));
    bufp->fullDouble(oldp+1293,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0xcU]));
    bufp->fullDouble(oldp+1295,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1297,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1299,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1301,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0xdU]));
    bufp->fullDouble(oldp+1303,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0xdU]));
    bufp->fullDouble(oldp+1305,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1307,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1309,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1311,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0xeU]));
    bufp->fullDouble(oldp+1313,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0xeU]));
    bufp->fullDouble(oldp+1315,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1317,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1319,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1321,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0xfU]));
    bufp->fullDouble(oldp+1323,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0xfU]));
    bufp->fullDouble(oldp+1325,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1327,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1329,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1331,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x10U]));
    bufp->fullDouble(oldp+1333,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x10U]));
    bufp->fullDouble(oldp+1335,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1337,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1339,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1341,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x11U]));
    bufp->fullDouble(oldp+1343,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x11U]));
    bufp->fullDouble(oldp+1345,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1347,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1349,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1351,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x12U]));
    bufp->fullDouble(oldp+1353,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x12U]));
    bufp->fullDouble(oldp+1355,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1357,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1359,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1361,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x13U]));
    bufp->fullDouble(oldp+1363,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x13U]));
    bufp->fullDouble(oldp+1365,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1367,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1369,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1371,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][1U]));
    bufp->fullDouble(oldp+1373,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][1U]));
    bufp->fullDouble(oldp+1375,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1377,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1379,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1381,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x14U]));
    bufp->fullDouble(oldp+1383,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x14U]));
    bufp->fullDouble(oldp+1385,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1387,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1389,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1391,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x15U]));
    bufp->fullDouble(oldp+1393,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x15U]));
    bufp->fullDouble(oldp+1395,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1397,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1399,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1401,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x16U]));
    bufp->fullDouble(oldp+1403,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x16U]));
    bufp->fullDouble(oldp+1405,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1407,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1409,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1411,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x17U]));
    bufp->fullDouble(oldp+1413,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x17U]));
    bufp->fullDouble(oldp+1415,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1417,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1419,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1421,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x18U]));
    bufp->fullDouble(oldp+1423,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x18U]));
    bufp->fullDouble(oldp+1425,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1427,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1429,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1431,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x19U]));
    bufp->fullDouble(oldp+1433,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x19U]));
    bufp->fullDouble(oldp+1435,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1437,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1439,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1441,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x1aU]));
    bufp->fullDouble(oldp+1443,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x1aU]));
    bufp->fullDouble(oldp+1445,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1447,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1449,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1451,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x1bU]));
    bufp->fullDouble(oldp+1453,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x1bU]));
    bufp->fullDouble(oldp+1455,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1457,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1459,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1461,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x1cU]));
    bufp->fullDouble(oldp+1463,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x1cU]));
    bufp->fullDouble(oldp+1465,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1467,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1469,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1471,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x1dU]));
    bufp->fullDouble(oldp+1473,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x1dU]));
    bufp->fullDouble(oldp+1475,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1477,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1479,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1481,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][2U]));
    bufp->fullDouble(oldp+1483,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][2U]));
    bufp->fullDouble(oldp+1485,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1487,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1489,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1491,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][0x1eU]));
    bufp->fullDouble(oldp+1493,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x1eU]));
    bufp->fullDouble(oldp+1495,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1497,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1499,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1501,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][3U]));
    bufp->fullDouble(oldp+1503,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][3U]));
    bufp->fullDouble(oldp+1505,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1507,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1509,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1511,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][4U]));
    bufp->fullDouble(oldp+1513,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][4U]));
    bufp->fullDouble(oldp+1515,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1517,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1519,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1521,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][5U]));
    bufp->fullDouble(oldp+1523,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][5U]));
    bufp->fullDouble(oldp+1525,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1527,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1529,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1531,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][6U]));
    bufp->fullDouble(oldp+1533,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][6U]));
    bufp->fullDouble(oldp+1535,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1537,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1539,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1541,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][7U]));
    bufp->fullDouble(oldp+1543,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][7U]));
    bufp->fullDouble(oldp+1545,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1547,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1549,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1551,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][8U]));
    bufp->fullDouble(oldp+1553,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][8U]));
    bufp->fullDouble(oldp+1555,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1557,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1559,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1561,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [0U][9U]));
    bufp->fullDouble(oldp+1563,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][9U]));
    bufp->fullDouble(oldp+1565,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1567,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1569,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1571,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [0U][0x1fU]));
    bufp->fullDouble(oldp+1573,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                                [1U]));
    bufp->fullDouble(oldp+1575,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_out));
    bufp->fullDouble(oldp+1577,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_sum_out));
    bufp->fullDouble(oldp+1579,((vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                 [0U][0x1fU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                                 [1U])));
    bufp->fullDouble(oldp+1581,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0U]));
    bufp->fullDouble(oldp+1583,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0U]));
    bufp->fullDouble(oldp+1585,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1587,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1589,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1591,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0xaU]));
    bufp->fullDouble(oldp+1593,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0xaU]));
    bufp->fullDouble(oldp+1595,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1597,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1599,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1601,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0xbU]));
    bufp->fullDouble(oldp+1603,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0xbU]));
    bufp->fullDouble(oldp+1605,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1607,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1609,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1611,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0xcU]));
    bufp->fullDouble(oldp+1613,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0xcU]));
    bufp->fullDouble(oldp+1615,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1617,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1619,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1621,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0xdU]));
    bufp->fullDouble(oldp+1623,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0xdU]));
    bufp->fullDouble(oldp+1625,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1627,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1629,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1631,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0xeU]));
    bufp->fullDouble(oldp+1633,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0xeU]));
    bufp->fullDouble(oldp+1635,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1637,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1639,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1641,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0xfU]));
    bufp->fullDouble(oldp+1643,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0xfU]));
    bufp->fullDouble(oldp+1645,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1647,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1649,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1651,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x10U]));
    bufp->fullDouble(oldp+1653,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x10U]));
    bufp->fullDouble(oldp+1655,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1657,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1659,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1661,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x11U]));
    bufp->fullDouble(oldp+1663,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x11U]));
    bufp->fullDouble(oldp+1665,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1667,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1669,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1671,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x12U]));
    bufp->fullDouble(oldp+1673,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x12U]));
    bufp->fullDouble(oldp+1675,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1677,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1679,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1681,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x13U]));
    bufp->fullDouble(oldp+1683,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x13U]));
    bufp->fullDouble(oldp+1685,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1687,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1689,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1691,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][1U]));
    bufp->fullDouble(oldp+1693,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][1U]));
    bufp->fullDouble(oldp+1695,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1697,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1699,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1701,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x14U]));
    bufp->fullDouble(oldp+1703,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x14U]));
    bufp->fullDouble(oldp+1705,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1707,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1709,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1711,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x15U]));
    bufp->fullDouble(oldp+1713,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x15U]));
    bufp->fullDouble(oldp+1715,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1717,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1719,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1721,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x16U]));
    bufp->fullDouble(oldp+1723,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x16U]));
    bufp->fullDouble(oldp+1725,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1727,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1729,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1731,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x17U]));
    bufp->fullDouble(oldp+1733,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x17U]));
    bufp->fullDouble(oldp+1735,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1737,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1739,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1741,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x18U]));
    bufp->fullDouble(oldp+1743,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x18U]));
    bufp->fullDouble(oldp+1745,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1747,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1749,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1751,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x19U]));
    bufp->fullDouble(oldp+1753,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x19U]));
    bufp->fullDouble(oldp+1755,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1757,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1759,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1761,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x1aU]));
    bufp->fullDouble(oldp+1763,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x1aU]));
    bufp->fullDouble(oldp+1765,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1767,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1769,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1771,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x1bU]));
    bufp->fullDouble(oldp+1773,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x1bU]));
    bufp->fullDouble(oldp+1775,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1777,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1779,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1781,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x1cU]));
    bufp->fullDouble(oldp+1783,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x1cU]));
    bufp->fullDouble(oldp+1785,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1787,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1789,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1791,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x1dU]));
    bufp->fullDouble(oldp+1793,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x1dU]));
    bufp->fullDouble(oldp+1795,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1797,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1799,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1801,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][2U]));
    bufp->fullDouble(oldp+1803,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][2U]));
    bufp->fullDouble(oldp+1805,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1807,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1809,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1811,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][0x1eU]));
    bufp->fullDouble(oldp+1813,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x1eU]));
    bufp->fullDouble(oldp+1815,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1817,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1819,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1821,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][3U]));
    bufp->fullDouble(oldp+1823,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][3U]));
    bufp->fullDouble(oldp+1825,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1827,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1829,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1831,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][4U]));
    bufp->fullDouble(oldp+1833,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][4U]));
    bufp->fullDouble(oldp+1835,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1837,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1839,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1841,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][5U]));
    bufp->fullDouble(oldp+1843,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][5U]));
    bufp->fullDouble(oldp+1845,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1847,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1849,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1851,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][6U]));
    bufp->fullDouble(oldp+1853,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][6U]));
    bufp->fullDouble(oldp+1855,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1857,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1859,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1861,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][7U]));
    bufp->fullDouble(oldp+1863,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][7U]));
    bufp->fullDouble(oldp+1865,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1867,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1869,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1871,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][8U]));
    bufp->fullDouble(oldp+1873,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][8U]));
    bufp->fullDouble(oldp+1875,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1877,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1879,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1881,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [1U][9U]));
    bufp->fullDouble(oldp+1883,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][9U]));
    bufp->fullDouble(oldp+1885,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1887,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1889,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1891,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [1U][0x1fU]));
    bufp->fullDouble(oldp+1893,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                                [2U]));
    bufp->fullDouble(oldp+1895,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_out));
    bufp->fullDouble(oldp+1897,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_sum_out));
    bufp->fullDouble(oldp+1899,((vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                 [1U][0x1fU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                                 [2U])));
    bufp->fullDouble(oldp+1901,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0U]));
    bufp->fullDouble(oldp+1903,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0U]));
    bufp->fullDouble(oldp+1905,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1907,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1909,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1911,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0xaU]));
    bufp->fullDouble(oldp+1913,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0xaU]));
    bufp->fullDouble(oldp+1915,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1917,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1919,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1921,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0xbU]));
    bufp->fullDouble(oldp+1923,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0xbU]));
    bufp->fullDouble(oldp+1925,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1927,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1929,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1931,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0xcU]));
    bufp->fullDouble(oldp+1933,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0xcU]));
    bufp->fullDouble(oldp+1935,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1937,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1939,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1941,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0xdU]));
    bufp->fullDouble(oldp+1943,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0xdU]));
    bufp->fullDouble(oldp+1945,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1947,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1949,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1951,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0xeU]));
    bufp->fullDouble(oldp+1953,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0xeU]));
    bufp->fullDouble(oldp+1955,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1957,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1959,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1961,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0xfU]));
    bufp->fullDouble(oldp+1963,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0xfU]));
    bufp->fullDouble(oldp+1965,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1967,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1969,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1971,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x10U]));
    bufp->fullDouble(oldp+1973,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x10U]));
    bufp->fullDouble(oldp+1975,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1977,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1979,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1981,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x11U]));
    bufp->fullDouble(oldp+1983,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x11U]));
    bufp->fullDouble(oldp+1985,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1987,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1989,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+1991,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x12U]));
    bufp->fullDouble(oldp+1993,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x12U]));
    bufp->fullDouble(oldp+1995,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+1997,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+1999,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2001,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x13U]));
    bufp->fullDouble(oldp+2003,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x13U]));
    bufp->fullDouble(oldp+2005,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2007,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2009,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2011,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][1U]));
    bufp->fullDouble(oldp+2013,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][1U]));
    bufp->fullDouble(oldp+2015,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2017,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2019,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2021,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x14U]));
    bufp->fullDouble(oldp+2023,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x14U]));
    bufp->fullDouble(oldp+2025,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2027,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2029,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2031,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x15U]));
    bufp->fullDouble(oldp+2033,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x15U]));
    bufp->fullDouble(oldp+2035,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2037,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2039,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2041,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x16U]));
    bufp->fullDouble(oldp+2043,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x16U]));
    bufp->fullDouble(oldp+2045,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2047,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2049,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2051,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x17U]));
    bufp->fullDouble(oldp+2053,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x17U]));
    bufp->fullDouble(oldp+2055,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2057,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2059,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2061,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x18U]));
    bufp->fullDouble(oldp+2063,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x18U]));
    bufp->fullDouble(oldp+2065,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2067,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2069,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2071,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x19U]));
    bufp->fullDouble(oldp+2073,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x19U]));
    bufp->fullDouble(oldp+2075,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2077,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2079,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2081,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x1aU]));
    bufp->fullDouble(oldp+2083,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x1aU]));
    bufp->fullDouble(oldp+2085,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2087,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2089,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2091,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x1bU]));
    bufp->fullDouble(oldp+2093,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x1bU]));
    bufp->fullDouble(oldp+2095,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2097,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2099,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2101,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x1cU]));
    bufp->fullDouble(oldp+2103,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x1cU]));
    bufp->fullDouble(oldp+2105,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2107,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2109,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2111,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x1dU]));
    bufp->fullDouble(oldp+2113,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x1dU]));
    bufp->fullDouble(oldp+2115,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2117,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2119,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2121,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][2U]));
    bufp->fullDouble(oldp+2123,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][2U]));
    bufp->fullDouble(oldp+2125,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2127,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2129,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2131,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][0x1eU]));
    bufp->fullDouble(oldp+2133,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x1eU]));
    bufp->fullDouble(oldp+2135,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2137,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2139,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2141,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][3U]));
    bufp->fullDouble(oldp+2143,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][3U]));
    bufp->fullDouble(oldp+2145,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2147,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2149,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2151,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][4U]));
    bufp->fullDouble(oldp+2153,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][4U]));
    bufp->fullDouble(oldp+2155,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2157,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2159,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2161,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][5U]));
    bufp->fullDouble(oldp+2163,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][5U]));
    bufp->fullDouble(oldp+2165,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2167,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2169,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2171,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][6U]));
    bufp->fullDouble(oldp+2173,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][6U]));
    bufp->fullDouble(oldp+2175,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2177,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2179,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2181,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][7U]));
    bufp->fullDouble(oldp+2183,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][7U]));
    bufp->fullDouble(oldp+2185,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2187,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2189,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2191,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][8U]));
    bufp->fullDouble(oldp+2193,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][8U]));
    bufp->fullDouble(oldp+2195,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2197,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2199,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2201,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [2U][9U]));
    bufp->fullDouble(oldp+2203,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][9U]));
    bufp->fullDouble(oldp+2205,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2207,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2209,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2211,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [2U][0x1fU]));
    bufp->fullDouble(oldp+2213,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                                [3U]));
    bufp->fullDouble(oldp+2215,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_out));
    bufp->fullDouble(oldp+2217,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_sum_out));
    bufp->fullDouble(oldp+2219,((vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                 [2U][0x1fU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                                 [3U])));
    bufp->fullDouble(oldp+2221,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0U]));
    bufp->fullDouble(oldp+2223,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0U]));
    bufp->fullDouble(oldp+2225,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2227,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2229,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2231,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0xaU]));
    bufp->fullDouble(oldp+2233,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0xaU]));
    bufp->fullDouble(oldp+2235,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2237,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2239,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2241,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0xbU]));
    bufp->fullDouble(oldp+2243,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0xbU]));
    bufp->fullDouble(oldp+2245,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2247,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2249,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2251,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0xcU]));
    bufp->fullDouble(oldp+2253,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0xcU]));
    bufp->fullDouble(oldp+2255,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2257,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2259,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2261,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0xdU]));
    bufp->fullDouble(oldp+2263,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0xdU]));
    bufp->fullDouble(oldp+2265,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2267,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2269,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2271,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0xeU]));
    bufp->fullDouble(oldp+2273,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0xeU]));
    bufp->fullDouble(oldp+2275,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2277,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2279,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2281,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0xfU]));
    bufp->fullDouble(oldp+2283,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0xfU]));
    bufp->fullDouble(oldp+2285,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2287,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2289,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2291,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x10U]));
    bufp->fullDouble(oldp+2293,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x10U]));
    bufp->fullDouble(oldp+2295,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2297,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2299,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2301,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x11U]));
    bufp->fullDouble(oldp+2303,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x11U]));
    bufp->fullDouble(oldp+2305,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2307,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2309,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2311,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x12U]));
    bufp->fullDouble(oldp+2313,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x12U]));
    bufp->fullDouble(oldp+2315,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2317,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2319,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2321,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x13U]));
    bufp->fullDouble(oldp+2323,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x13U]));
    bufp->fullDouble(oldp+2325,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2327,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2329,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2331,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+2333,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+2335,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2337,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2339,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2341,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x14U]));
    bufp->fullDouble(oldp+2343,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x14U]));
    bufp->fullDouble(oldp+2345,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2347,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2349,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2351,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x15U]));
    bufp->fullDouble(oldp+2353,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x15U]));
    bufp->fullDouble(oldp+2355,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2357,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2359,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2361,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x16U]));
    bufp->fullDouble(oldp+2363,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x16U]));
    bufp->fullDouble(oldp+2365,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2367,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2369,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2371,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x17U]));
    bufp->fullDouble(oldp+2373,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x17U]));
    bufp->fullDouble(oldp+2375,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2377,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2379,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2381,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x18U]));
    bufp->fullDouble(oldp+2383,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x18U]));
    bufp->fullDouble(oldp+2385,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2387,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2389,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2391,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x19U]));
    bufp->fullDouble(oldp+2393,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x19U]));
    bufp->fullDouble(oldp+2395,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2397,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2399,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2401,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x1aU]));
    bufp->fullDouble(oldp+2403,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x1aU]));
    bufp->fullDouble(oldp+2405,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2407,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2409,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2411,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x1bU]));
    bufp->fullDouble(oldp+2413,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x1bU]));
    bufp->fullDouble(oldp+2415,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2417,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2419,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2421,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x1cU]));
    bufp->fullDouble(oldp+2423,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x1cU]));
    bufp->fullDouble(oldp+2425,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2427,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2429,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2431,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x1dU]));
    bufp->fullDouble(oldp+2433,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x1dU]));
    bufp->fullDouble(oldp+2435,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2437,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2439,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2441,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+2443,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+2445,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2447,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2449,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2451,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][0x1eU]));
    bufp->fullDouble(oldp+2453,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x1eU]));
    bufp->fullDouble(oldp+2455,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2457,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2459,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2461,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+2463,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+2465,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2467,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2469,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2471,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][4U]));
    bufp->fullDouble(oldp+2473,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][4U]));
    bufp->fullDouble(oldp+2475,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2477,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2479,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2481,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][5U]));
    bufp->fullDouble(oldp+2483,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][5U]));
    bufp->fullDouble(oldp+2485,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2487,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2489,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2491,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][6U]));
    bufp->fullDouble(oldp+2493,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][6U]));
    bufp->fullDouble(oldp+2495,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2497,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2499,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2501,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][7U]));
    bufp->fullDouble(oldp+2503,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][7U]));
    bufp->fullDouble(oldp+2505,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2507,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2509,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2511,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][8U]));
    bufp->fullDouble(oldp+2513,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][8U]));
    bufp->fullDouble(oldp+2515,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2517,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2519,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2521,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_pow_wire
                                [3U][9U]));
    bufp->fullDouble(oldp+2523,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][9U]));
    bufp->fullDouble(oldp+2525,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__x_out));
    bufp->fullDouble(oldp+2527,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out));
    bufp->fullDouble(oldp+2529,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out));
    bufp->fullDouble(oldp+2531,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                [3U][0x1fU]));
    bufp->fullDouble(oldp+2533,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                                [4U]));
    bufp->fullDouble(oldp+2535,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_out));
    bufp->fullDouble(oldp+2537,(vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_sum_out));
    bufp->fullDouble(oldp+2539,((vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                                 [3U][0x1fU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                                 [4U])));
    bufp->fullDouble(oldp+2541,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[0]));
    bufp->fullDouble(oldp+2543,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[1]));
    bufp->fullDouble(oldp+2545,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[2]));
    bufp->fullDouble(oldp+2547,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[3]));
    bufp->fullDouble(oldp+2549,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire[4]));
    bufp->fullDouble(oldp+2551,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__dummy_pe__data_out));
    bufp->fullDouble(oldp+2553,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                                [0U]));
    bufp->fullDouble(oldp+2555,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                                [1U]));
    bufp->fullDouble(oldp+2557,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__x_normed_out));
    bufp->fullDouble(oldp+2559,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__norm_out));
    bufp->fullDouble(oldp+2561,((vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                                 [0U] / vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                                 [1U])));
    bufp->fullDouble(oldp+2563,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                                [1U]));
    bufp->fullDouble(oldp+2565,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                                [2U]));
    bufp->fullDouble(oldp+2567,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__x_normed_out));
    bufp->fullDouble(oldp+2569,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__norm_out));
    bufp->fullDouble(oldp+2571,((vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                                 [1U] / vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                                 [2U])));
    bufp->fullDouble(oldp+2573,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                                [2U]));
    bufp->fullDouble(oldp+2575,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                                [3U]));
    bufp->fullDouble(oldp+2577,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__x_normed_out));
    bufp->fullDouble(oldp+2579,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__norm_out));
    bufp->fullDouble(oldp+2581,((vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                                 [2U] / vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                                 [3U])));
    bufp->fullDouble(oldp+2583,(vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                                [3U]));
    bufp->fullDouble(oldp+2585,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                                [4U]));
    bufp->fullDouble(oldp+2587,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__x_normed_out));
    bufp->fullDouble(oldp+2589,(vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__norm_out));
    bufp->fullDouble(oldp+2591,((vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                                 [3U] / vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                                 [4U])));
    bufp->fullIData(oldp+2593,(vlSelfRef.attentionSystolicArray__DOT__unnamedblk1__DOT__j),32);
    bufp->fullIData(oldp+2594,(vlSelfRef.attentionSystolicArray__DOT__unnamedblk11__DOT__i),32);
    bufp->fullIData(oldp+2595,(vlSelfRef.attentionSystolicArray__DOT__unnamedblk11__DOT__unnamedblk12__DOT__j),32);
    bufp->fullIData(oldp+2596,(vlSelfRef.attentionSystolicArray__DOT__unnamedblk13__DOT__i),32);
    bufp->fullIData(oldp+2597,(vlSelfRef.attentionSystolicArray__DOT__unnamedblk13__DOT__unnamedblk14__DOT__j),32);
    bufp->fullIData(oldp+2598,(vlSelfRef.attentionSystolicArray__DOT__unnamedblk15__DOT__i),32);
    bufp->fullIData(oldp+2599,(vlSelfRef.attentionSystolicArray__DOT__unnamedblk15__DOT__unnamedblk16__DOT__j),32);
    bufp->fullIData(oldp+2600,(vlSelfRef.attentionSystolicArray__DOT__unnamedblk2__DOT__j),32);
    bufp->fullIData(oldp+2601,(vlSelfRef.attentionSystolicArray__DOT__unnamedblk9__DOT__i),32);
    bufp->fullIData(oldp+2602,(vlSelfRef.attentionSystolicArray__DOT__unnamedblk9__DOT__unnamedblk10__DOT__j),32);
    bufp->fullDouble(oldp+2603,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.data_in[0]));
    bufp->fullDouble(oldp+2605,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.data_in[1]));
    bufp->fullDouble(oldp+2607,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.data_in[2]));
    bufp->fullDouble(oldp+2609,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.data_in[3]));
    bufp->fullDouble(oldp+2611,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.result_out[0]));
    bufp->fullDouble(oldp+2613,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.result_out[1]));
    bufp->fullDouble(oldp+2615,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.result_out[2]));
    bufp->fullDouble(oldp+2617,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.result_out[3]));
    bufp->fullDouble(oldp+2619,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [0U][0U]));
    bufp->fullDouble(oldp+2621,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [0U][1U]));
    bufp->fullDouble(oldp+2623,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [0U][2U]));
    bufp->fullDouble(oldp+2625,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [0U][3U]));
    bufp->fullDouble(oldp+2627,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [1U][0U]));
    bufp->fullDouble(oldp+2629,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [1U][1U]));
    bufp->fullDouble(oldp+2631,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [1U][2U]));
    bufp->fullDouble(oldp+2633,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [1U][3U]));
    bufp->fullDouble(oldp+2635,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [2U][0U]));
    bufp->fullDouble(oldp+2637,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [2U][1U]));
    bufp->fullDouble(oldp+2639,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [2U][2U]));
    bufp->fullDouble(oldp+2641,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [2U][3U]));
    bufp->fullDouble(oldp+2643,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [3U][0U]));
    bufp->fullDouble(oldp+2645,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+2647,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+2649,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+2651,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [4U][0U]));
    bufp->fullDouble(oldp+2653,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [4U][1U]));
    bufp->fullDouble(oldp+2655,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [4U][2U]));
    bufp->fullDouble(oldp+2657,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [4U][3U]));
    bufp->fullDouble(oldp+2659,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [0U][0U]));
    bufp->fullDouble(oldp+2661,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [0U][1U]));
    bufp->fullDouble(oldp+2663,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [0U][2U]));
    bufp->fullDouble(oldp+2665,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [0U][3U]));
    bufp->fullDouble(oldp+2667,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [0U][4U]));
    bufp->fullDouble(oldp+2669,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [1U][0U]));
    bufp->fullDouble(oldp+2671,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [1U][1U]));
    bufp->fullDouble(oldp+2673,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [1U][2U]));
    bufp->fullDouble(oldp+2675,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [1U][3U]));
    bufp->fullDouble(oldp+2677,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [1U][4U]));
    bufp->fullDouble(oldp+2679,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [2U][0U]));
    bufp->fullDouble(oldp+2681,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [2U][1U]));
    bufp->fullDouble(oldp+2683,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [2U][2U]));
    bufp->fullDouble(oldp+2685,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [2U][3U]));
    bufp->fullDouble(oldp+2687,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [2U][4U]));
    bufp->fullDouble(oldp+2689,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [3U][0U]));
    bufp->fullDouble(oldp+2691,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+2693,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+2695,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+2697,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [3U][4U]));
    bufp->fullDouble(oldp+2699,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [1U][0U]));
    bufp->fullDouble(oldp+2701,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [0U][0U]));
    bufp->fullDouble(oldp+2703,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2705,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2707,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2709,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2711,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2713,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [1U][1U]));
    bufp->fullDouble(oldp+2715,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [0U][1U]));
    bufp->fullDouble(oldp+2717,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2719,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2721,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2723,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2725,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2727,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [1U][2U]));
    bufp->fullDouble(oldp+2729,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [0U][2U]));
    bufp->fullDouble(oldp+2731,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2733,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2735,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2737,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2739,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2741,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [1U][3U]));
    bufp->fullDouble(oldp+2743,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [0U][3U]));
    bufp->fullDouble(oldp+2745,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2747,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2749,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2751,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2753,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2755,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [2U][0U]));
    bufp->fullDouble(oldp+2757,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [1U][0U]));
    bufp->fullDouble(oldp+2759,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2761,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2763,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2765,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2767,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2769,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [2U][1U]));
    bufp->fullDouble(oldp+2771,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [1U][1U]));
    bufp->fullDouble(oldp+2773,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2775,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2777,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2779,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2781,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2783,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [2U][2U]));
    bufp->fullDouble(oldp+2785,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [1U][2U]));
    bufp->fullDouble(oldp+2787,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2789,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2791,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2793,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2795,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2797,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [2U][3U]));
    bufp->fullDouble(oldp+2799,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [1U][3U]));
    bufp->fullDouble(oldp+2801,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2803,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2805,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2807,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2809,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2811,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [3U][0U]));
    bufp->fullDouble(oldp+2813,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [2U][0U]));
    bufp->fullDouble(oldp+2815,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2817,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2819,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2821,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2823,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2825,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+2827,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [2U][1U]));
    bufp->fullDouble(oldp+2829,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2831,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2833,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2835,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2837,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2839,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+2841,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [2U][2U]));
    bufp->fullDouble(oldp+2843,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2845,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2847,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2849,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2851,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2853,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+2855,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [2U][3U]));
    bufp->fullDouble(oldp+2857,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2859,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2861,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2863,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2865,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2867,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [4U][0U]));
    bufp->fullDouble(oldp+2869,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [3U][0U]));
    bufp->fullDouble(oldp+2871,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2873,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2875,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2877,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2879,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2881,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [4U][1U]));
    bufp->fullDouble(oldp+2883,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+2885,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2887,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2889,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2891,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2893,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2895,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [4U][2U]));
    bufp->fullDouble(oldp+2897,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+2899,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2901,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2903,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2905,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2907,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2909,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__data_wire
                                [4U][3U]));
    bufp->fullDouble(oldp+2911,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__partial_sum_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+2913,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+2915,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+2917,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+2919,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+2921,(vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+2923,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.data_in[0]));
    bufp->fullDouble(oldp+2925,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.data_in[1]));
    bufp->fullDouble(oldp+2927,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.data_in[2]));
    bufp->fullDouble(oldp+2929,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.data_in[3]));
    bufp->fullDouble(oldp+2931,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.result_out[0]));
    bufp->fullDouble(oldp+2933,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.result_out[1]));
    bufp->fullDouble(oldp+2935,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.result_out[2]));
    bufp->fullDouble(oldp+2937,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.result_out[3]));
    bufp->fullDouble(oldp+2939,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [0U][0U]));
    bufp->fullDouble(oldp+2941,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [0U][1U]));
    bufp->fullDouble(oldp+2943,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [0U][2U]));
    bufp->fullDouble(oldp+2945,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [0U][3U]));
    bufp->fullDouble(oldp+2947,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [1U][0U]));
    bufp->fullDouble(oldp+2949,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [1U][1U]));
    bufp->fullDouble(oldp+2951,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [1U][2U]));
    bufp->fullDouble(oldp+2953,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [1U][3U]));
    bufp->fullDouble(oldp+2955,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [2U][0U]));
    bufp->fullDouble(oldp+2957,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [2U][1U]));
    bufp->fullDouble(oldp+2959,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [2U][2U]));
    bufp->fullDouble(oldp+2961,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [2U][3U]));
    bufp->fullDouble(oldp+2963,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [3U][0U]));
    bufp->fullDouble(oldp+2965,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+2967,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+2969,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+2971,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [4U][0U]));
    bufp->fullDouble(oldp+2973,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [4U][1U]));
    bufp->fullDouble(oldp+2975,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [4U][2U]));
    bufp->fullDouble(oldp+2977,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [4U][3U]));
    bufp->fullDouble(oldp+2979,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [0U][0U]));
    bufp->fullDouble(oldp+2981,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [0U][1U]));
    bufp->fullDouble(oldp+2983,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [0U][2U]));
    bufp->fullDouble(oldp+2985,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [0U][3U]));
    bufp->fullDouble(oldp+2987,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [0U][4U]));
    bufp->fullDouble(oldp+2989,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [1U][0U]));
    bufp->fullDouble(oldp+2991,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [1U][1U]));
    bufp->fullDouble(oldp+2993,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [1U][2U]));
    bufp->fullDouble(oldp+2995,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [1U][3U]));
    bufp->fullDouble(oldp+2997,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [1U][4U]));
    bufp->fullDouble(oldp+2999,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [2U][0U]));
    bufp->fullDouble(oldp+3001,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [2U][1U]));
    bufp->fullDouble(oldp+3003,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [2U][2U]));
    bufp->fullDouble(oldp+3005,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [2U][3U]));
    bufp->fullDouble(oldp+3007,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [2U][4U]));
    bufp->fullDouble(oldp+3009,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [3U][0U]));
    bufp->fullDouble(oldp+3011,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+3013,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+3015,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+3017,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [3U][4U]));
    bufp->fullDouble(oldp+3019,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [1U][0U]));
    bufp->fullDouble(oldp+3021,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [0U][0U]));
    bufp->fullDouble(oldp+3023,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3025,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3027,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3029,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3031,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3033,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [1U][1U]));
    bufp->fullDouble(oldp+3035,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [0U][1U]));
    bufp->fullDouble(oldp+3037,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3039,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3041,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3043,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3045,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3047,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [1U][2U]));
    bufp->fullDouble(oldp+3049,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [0U][2U]));
    bufp->fullDouble(oldp+3051,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3053,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3055,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3057,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3059,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3061,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [1U][3U]));
    bufp->fullDouble(oldp+3063,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [0U][3U]));
    bufp->fullDouble(oldp+3065,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3067,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3069,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3071,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3073,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3075,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [2U][0U]));
    bufp->fullDouble(oldp+3077,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [1U][0U]));
    bufp->fullDouble(oldp+3079,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3081,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3083,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3085,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3087,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3089,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [2U][1U]));
    bufp->fullDouble(oldp+3091,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [1U][1U]));
    bufp->fullDouble(oldp+3093,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3095,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3097,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3099,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3101,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3103,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [2U][2U]));
    bufp->fullDouble(oldp+3105,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [1U][2U]));
    bufp->fullDouble(oldp+3107,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3109,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3111,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3113,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3115,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3117,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [2U][3U]));
    bufp->fullDouble(oldp+3119,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [1U][3U]));
    bufp->fullDouble(oldp+3121,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3123,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3125,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3127,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3129,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3131,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [3U][0U]));
    bufp->fullDouble(oldp+3133,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [2U][0U]));
    bufp->fullDouble(oldp+3135,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3137,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3139,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3141,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3143,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3145,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+3147,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [2U][1U]));
    bufp->fullDouble(oldp+3149,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3151,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3153,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3155,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3157,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3159,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+3161,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [2U][2U]));
    bufp->fullDouble(oldp+3163,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3165,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3167,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3169,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3171,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3173,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+3175,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [2U][3U]));
    bufp->fullDouble(oldp+3177,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3179,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3181,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3183,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3185,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3187,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [4U][0U]));
    bufp->fullDouble(oldp+3189,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [3U][0U]));
    bufp->fullDouble(oldp+3191,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3193,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3195,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3197,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3199,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3201,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [4U][1U]));
    bufp->fullDouble(oldp+3203,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [3U][1U]));
    bufp->fullDouble(oldp+3205,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3207,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3209,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3211,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3213,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3215,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [4U][2U]));
    bufp->fullDouble(oldp+3217,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [3U][2U]));
    bufp->fullDouble(oldp+3219,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3221,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3223,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3225,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3227,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__add_result));
    bufp->fullDouble(oldp+3229,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__data_wire
                                [4U][3U]));
    bufp->fullDouble(oldp+3231,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__partial_sum_wire
                                [3U][3U]));
    bufp->fullDouble(oldp+3233,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__data_out));
    bufp->fullDouble(oldp+3235,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__partial_sum_out));
    bufp->fullDouble(oldp+3237,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__weight));
    bufp->fullDouble(oldp+3239,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__mul_result));
    bufp->fullDouble(oldp+3241,(vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray.__PVT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__add_result));
    bufp->fullBit(oldp+3243,(vlSelfRef.clk));
    bufp->fullBit(oldp+3244,(vlSelfRef.reset));
    bufp->fullBit(oldp+3245,(vlSelfRef.valid_input));
    bufp->fullDouble(oldp+3246,(vlSelfRef.K_matrix[0U]
                                [0U]));
    bufp->fullDouble(oldp+3248,(vlSelfRef.K_matrix[0U]
                                [1U]));
    bufp->fullDouble(oldp+3250,(vlSelfRef.K_matrix[0U]
                                [2U]));
    bufp->fullDouble(oldp+3252,(vlSelfRef.K_matrix[0U]
                                [3U]));
    bufp->fullDouble(oldp+3254,(vlSelfRef.K_matrix[1U]
                                [0U]));
    bufp->fullDouble(oldp+3256,(vlSelfRef.K_matrix[1U]
                                [1U]));
    bufp->fullDouble(oldp+3258,(vlSelfRef.K_matrix[1U]
                                [2U]));
    bufp->fullDouble(oldp+3260,(vlSelfRef.K_matrix[1U]
                                [3U]));
    bufp->fullDouble(oldp+3262,(vlSelfRef.K_matrix[2U]
                                [0U]));
    bufp->fullDouble(oldp+3264,(vlSelfRef.K_matrix[2U]
                                [1U]));
    bufp->fullDouble(oldp+3266,(vlSelfRef.K_matrix[2U]
                                [2U]));
    bufp->fullDouble(oldp+3268,(vlSelfRef.K_matrix[2U]
                                [3U]));
    bufp->fullDouble(oldp+3270,(vlSelfRef.K_matrix[3U]
                                [0U]));
    bufp->fullDouble(oldp+3272,(vlSelfRef.K_matrix[3U]
                                [1U]));
    bufp->fullDouble(oldp+3274,(vlSelfRef.K_matrix[3U]
                                [2U]));
    bufp->fullDouble(oldp+3276,(vlSelfRef.K_matrix[3U]
                                [3U]));
    bufp->fullDouble(oldp+3278,(vlSelfRef.Q_matrix[0U]
                                [0U]));
    bufp->fullDouble(oldp+3280,(vlSelfRef.Q_matrix[0U]
                                [1U]));
    bufp->fullDouble(oldp+3282,(vlSelfRef.Q_matrix[0U]
                                [2U]));
    bufp->fullDouble(oldp+3284,(vlSelfRef.Q_matrix[0U]
                                [3U]));
    bufp->fullDouble(oldp+3286,(vlSelfRef.Q_matrix[1U]
                                [0U]));
    bufp->fullDouble(oldp+3288,(vlSelfRef.Q_matrix[1U]
                                [1U]));
    bufp->fullDouble(oldp+3290,(vlSelfRef.Q_matrix[1U]
                                [2U]));
    bufp->fullDouble(oldp+3292,(vlSelfRef.Q_matrix[1U]
                                [3U]));
    bufp->fullDouble(oldp+3294,(vlSelfRef.Q_matrix[2U]
                                [0U]));
    bufp->fullDouble(oldp+3296,(vlSelfRef.Q_matrix[2U]
                                [1U]));
    bufp->fullDouble(oldp+3298,(vlSelfRef.Q_matrix[2U]
                                [2U]));
    bufp->fullDouble(oldp+3300,(vlSelfRef.Q_matrix[2U]
                                [3U]));
    bufp->fullDouble(oldp+3302,(vlSelfRef.Q_matrix[3U]
                                [0U]));
    bufp->fullDouble(oldp+3304,(vlSelfRef.Q_matrix[3U]
                                [1U]));
    bufp->fullDouble(oldp+3306,(vlSelfRef.Q_matrix[3U]
                                [2U]));
    bufp->fullDouble(oldp+3308,(vlSelfRef.Q_matrix[3U]
                                [3U]));
    bufp->fullDouble(oldp+3310,(vlSelfRef.V_matrix[0U]
                                [0U]));
    bufp->fullDouble(oldp+3312,(vlSelfRef.V_matrix[0U]
                                [1U]));
    bufp->fullDouble(oldp+3314,(vlSelfRef.V_matrix[0U]
                                [2U]));
    bufp->fullDouble(oldp+3316,(vlSelfRef.V_matrix[0U]
                                [3U]));
    bufp->fullDouble(oldp+3318,(vlSelfRef.V_matrix[1U]
                                [0U]));
    bufp->fullDouble(oldp+3320,(vlSelfRef.V_matrix[1U]
                                [1U]));
    bufp->fullDouble(oldp+3322,(vlSelfRef.V_matrix[1U]
                                [2U]));
    bufp->fullDouble(oldp+3324,(vlSelfRef.V_matrix[1U]
                                [3U]));
    bufp->fullDouble(oldp+3326,(vlSelfRef.V_matrix[2U]
                                [0U]));
    bufp->fullDouble(oldp+3328,(vlSelfRef.V_matrix[2U]
                                [1U]));
    bufp->fullDouble(oldp+3330,(vlSelfRef.V_matrix[2U]
                                [2U]));
    bufp->fullDouble(oldp+3332,(vlSelfRef.V_matrix[2U]
                                [3U]));
    bufp->fullDouble(oldp+3334,(vlSelfRef.V_matrix[3U]
                                [0U]));
    bufp->fullDouble(oldp+3336,(vlSelfRef.V_matrix[3U]
                                [1U]));
    bufp->fullDouble(oldp+3338,(vlSelfRef.V_matrix[3U]
                                [2U]));
    bufp->fullDouble(oldp+3340,(vlSelfRef.V_matrix[3U]
                                [3U]));
    bufp->fullDouble(oldp+3342,(vlSelfRef.softmax_temp));
    bufp->fullDouble(oldp+3344,(vlSelfRef.Q_mult_K[0U]
                                [0U]));
    bufp->fullDouble(oldp+3346,(vlSelfRef.Q_mult_K[0U]
                                [1U]));
    bufp->fullDouble(oldp+3348,(vlSelfRef.Q_mult_K[0U]
                                [2U]));
    bufp->fullDouble(oldp+3350,(vlSelfRef.Q_mult_K[0U]
                                [3U]));
    bufp->fullDouble(oldp+3352,(vlSelfRef.Q_mult_K[1U]
                                [0U]));
    bufp->fullDouble(oldp+3354,(vlSelfRef.Q_mult_K[1U]
                                [1U]));
    bufp->fullDouble(oldp+3356,(vlSelfRef.Q_mult_K[1U]
                                [2U]));
    bufp->fullDouble(oldp+3358,(vlSelfRef.Q_mult_K[1U]
                                [3U]));
    bufp->fullDouble(oldp+3360,(vlSelfRef.Q_mult_K[2U]
                                [0U]));
    bufp->fullDouble(oldp+3362,(vlSelfRef.Q_mult_K[2U]
                                [1U]));
    bufp->fullDouble(oldp+3364,(vlSelfRef.Q_mult_K[2U]
                                [2U]));
    bufp->fullDouble(oldp+3366,(vlSelfRef.Q_mult_K[2U]
                                [3U]));
    bufp->fullDouble(oldp+3368,(vlSelfRef.Q_mult_K[3U]
                                [0U]));
    bufp->fullDouble(oldp+3370,(vlSelfRef.Q_mult_K[3U]
                                [1U]));
    bufp->fullDouble(oldp+3372,(vlSelfRef.Q_mult_K[3U]
                                [2U]));
    bufp->fullDouble(oldp+3374,(vlSelfRef.Q_mult_K[3U]
                                [3U]));
    bufp->fullDouble(oldp+3376,(vlSelfRef.exp_Q_mult_K
                                [0U][0U]));
    bufp->fullDouble(oldp+3378,(vlSelfRef.exp_Q_mult_K
                                [0U][1U]));
    bufp->fullDouble(oldp+3380,(vlSelfRef.exp_Q_mult_K
                                [0U][2U]));
    bufp->fullDouble(oldp+3382,(vlSelfRef.exp_Q_mult_K
                                [0U][3U]));
    bufp->fullDouble(oldp+3384,(vlSelfRef.exp_Q_mult_K
                                [1U][0U]));
    bufp->fullDouble(oldp+3386,(vlSelfRef.exp_Q_mult_K
                                [1U][1U]));
    bufp->fullDouble(oldp+3388,(vlSelfRef.exp_Q_mult_K
                                [1U][2U]));
    bufp->fullDouble(oldp+3390,(vlSelfRef.exp_Q_mult_K
                                [1U][3U]));
    bufp->fullDouble(oldp+3392,(vlSelfRef.exp_Q_mult_K
                                [2U][0U]));
    bufp->fullDouble(oldp+3394,(vlSelfRef.exp_Q_mult_K
                                [2U][1U]));
    bufp->fullDouble(oldp+3396,(vlSelfRef.exp_Q_mult_K
                                [2U][2U]));
    bufp->fullDouble(oldp+3398,(vlSelfRef.exp_Q_mult_K
                                [2U][3U]));
    bufp->fullDouble(oldp+3400,(vlSelfRef.exp_Q_mult_K
                                [3U][0U]));
    bufp->fullDouble(oldp+3402,(vlSelfRef.exp_Q_mult_K
                                [3U][1U]));
    bufp->fullDouble(oldp+3404,(vlSelfRef.exp_Q_mult_K
                                [3U][2U]));
    bufp->fullDouble(oldp+3406,(vlSelfRef.exp_Q_mult_K
                                [3U][3U]));
    bufp->fullDouble(oldp+3408,(vlSelfRef.exp_K_mult_Q_mult_V
                                [0U][0U]));
    bufp->fullDouble(oldp+3410,(vlSelfRef.exp_K_mult_Q_mult_V
                                [0U][1U]));
    bufp->fullDouble(oldp+3412,(vlSelfRef.exp_K_mult_Q_mult_V
                                [0U][2U]));
    bufp->fullDouble(oldp+3414,(vlSelfRef.exp_K_mult_Q_mult_V
                                [0U][3U]));
    bufp->fullDouble(oldp+3416,(vlSelfRef.exp_K_mult_Q_mult_V
                                [1U][0U]));
    bufp->fullDouble(oldp+3418,(vlSelfRef.exp_K_mult_Q_mult_V
                                [1U][1U]));
    bufp->fullDouble(oldp+3420,(vlSelfRef.exp_K_mult_Q_mult_V
                                [1U][2U]));
    bufp->fullDouble(oldp+3422,(vlSelfRef.exp_K_mult_Q_mult_V
                                [1U][3U]));
    bufp->fullDouble(oldp+3424,(vlSelfRef.exp_K_mult_Q_mult_V
                                [2U][0U]));
    bufp->fullDouble(oldp+3426,(vlSelfRef.exp_K_mult_Q_mult_V
                                [2U][1U]));
    bufp->fullDouble(oldp+3428,(vlSelfRef.exp_K_mult_Q_mult_V
                                [2U][2U]));
    bufp->fullDouble(oldp+3430,(vlSelfRef.exp_K_mult_Q_mult_V
                                [2U][3U]));
    bufp->fullDouble(oldp+3432,(vlSelfRef.exp_K_mult_Q_mult_V
                                [3U][0U]));
    bufp->fullDouble(oldp+3434,(vlSelfRef.exp_K_mult_Q_mult_V
                                [3U][1U]));
    bufp->fullDouble(oldp+3436,(vlSelfRef.exp_K_mult_Q_mult_V
                                [3U][2U]));
    bufp->fullDouble(oldp+3438,(vlSelfRef.exp_K_mult_Q_mult_V
                                [3U][3U]));
    bufp->fullDouble(oldp+3440,(vlSelfRef.attention
                                [0U][0U]));
    bufp->fullDouble(oldp+3442,(vlSelfRef.attention
                                [0U][1U]));
    bufp->fullDouble(oldp+3444,(vlSelfRef.attention
                                [0U][2U]));
    bufp->fullDouble(oldp+3446,(vlSelfRef.attention
                                [0U][3U]));
    bufp->fullDouble(oldp+3448,(vlSelfRef.attention
                                [1U][0U]));
    bufp->fullDouble(oldp+3450,(vlSelfRef.attention
                                [1U][1U]));
    bufp->fullDouble(oldp+3452,(vlSelfRef.attention
                                [1U][2U]));
    bufp->fullDouble(oldp+3454,(vlSelfRef.attention
                                [1U][3U]));
    bufp->fullDouble(oldp+3456,(vlSelfRef.attention
                                [2U][0U]));
    bufp->fullDouble(oldp+3458,(vlSelfRef.attention
                                [2U][1U]));
    bufp->fullDouble(oldp+3460,(vlSelfRef.attention
                                [2U][2U]));
    bufp->fullDouble(oldp+3462,(vlSelfRef.attention
                                [2U][3U]));
    bufp->fullDouble(oldp+3464,(vlSelfRef.attention
                                [3U][0U]));
    bufp->fullDouble(oldp+3466,(vlSelfRef.attention
                                [3U][1U]));
    bufp->fullDouble(oldp+3468,(vlSelfRef.attention
                                [3U][2U]));
    bufp->fullDouble(oldp+3470,(vlSelfRef.attention
                                [3U][3U]));
    bufp->fullBit(oldp+3472,(vlSelfRef.valid_result));
    bufp->fullDouble(oldp+3473,(vlSelfRef.factorial_arr
                                [1U]));
    bufp->fullDouble(oldp+3475,(vlSelfRef.factorial_arr
                                [0xbU]));
    bufp->fullDouble(oldp+3477,(vlSelfRef.factorial_arr
                                [0xcU]));
    bufp->fullDouble(oldp+3479,(vlSelfRef.factorial_arr
                                [0xdU]));
    bufp->fullDouble(oldp+3481,(vlSelfRef.factorial_arr
                                [0xeU]));
    bufp->fullDouble(oldp+3483,(vlSelfRef.factorial_arr
                                [0xfU]));
    bufp->fullDouble(oldp+3485,(vlSelfRef.factorial_arr
                                [0x10U]));
    bufp->fullDouble(oldp+3487,(vlSelfRef.factorial_arr
                                [0x11U]));
    bufp->fullDouble(oldp+3489,(vlSelfRef.factorial_arr
                                [0x12U]));
    bufp->fullDouble(oldp+3491,(vlSelfRef.factorial_arr
                                [0x13U]));
    bufp->fullDouble(oldp+3493,(vlSelfRef.factorial_arr
                                [0x14U]));
    bufp->fullDouble(oldp+3495,(vlSelfRef.factorial_arr
                                [2U]));
    bufp->fullDouble(oldp+3497,(vlSelfRef.factorial_arr
                                [0x15U]));
    bufp->fullDouble(oldp+3499,(vlSelfRef.factorial_arr
                                [0x16U]));
    bufp->fullDouble(oldp+3501,(vlSelfRef.factorial_arr
                                [0x17U]));
    bufp->fullDouble(oldp+3503,(vlSelfRef.factorial_arr
                                [0x18U]));
    bufp->fullDouble(oldp+3505,(vlSelfRef.factorial_arr
                                [0x19U]));
    bufp->fullDouble(oldp+3507,(vlSelfRef.factorial_arr
                                [0x1aU]));
    bufp->fullDouble(oldp+3509,(vlSelfRef.factorial_arr
                                [0x1bU]));
    bufp->fullDouble(oldp+3511,(vlSelfRef.factorial_arr
                                [0x1cU]));
    bufp->fullDouble(oldp+3513,(vlSelfRef.factorial_arr
                                [0x1dU]));
    bufp->fullDouble(oldp+3515,(vlSelfRef.factorial_arr
                                [0x1eU]));
    bufp->fullDouble(oldp+3517,(vlSelfRef.factorial_arr
                                [3U]));
    bufp->fullDouble(oldp+3519,(vlSelfRef.factorial_arr
                                [0x1fU]));
    bufp->fullDouble(oldp+3521,(vlSelfRef.factorial_arr
                                [4U]));
    bufp->fullDouble(oldp+3523,(vlSelfRef.factorial_arr
                                [5U]));
    bufp->fullDouble(oldp+3525,(vlSelfRef.factorial_arr
                                [6U]));
    bufp->fullDouble(oldp+3527,(vlSelfRef.factorial_arr
                                [7U]));
    bufp->fullDouble(oldp+3529,(vlSelfRef.factorial_arr
                                [8U]));
    bufp->fullDouble(oldp+3531,(vlSelfRef.factorial_arr
                                [9U]));
    bufp->fullDouble(oldp+3533,(vlSelfRef.factorial_arr
                                [0xaU]));
}
