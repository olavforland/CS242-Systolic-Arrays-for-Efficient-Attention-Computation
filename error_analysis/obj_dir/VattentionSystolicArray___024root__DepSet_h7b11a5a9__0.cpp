// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Design implementation internals
// See VattentionSystolicArray.h for the primary calling header

#include "VattentionSystolicArray__pch.h"
#include "VattentionSystolicArray__Syms.h"
#include "VattentionSystolicArray___024root.h"

#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__ico(VattentionSystolicArray___024root* vlSelf);
#endif  // VL_DEBUG

void VattentionSystolicArray___024root___eval_triggers__ico(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_triggers__ico\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelfRef.__VicoTriggered.set(0U, (IData)(vlSelfRef.__VicoFirstIteration));
#ifdef VL_DEBUG
    if (VL_UNLIKELY(vlSymsp->_vm_contextp__->debug())) {
        VattentionSystolicArray___024root___dump_triggers__ico(vlSelf);
    }
#endif
}

VL_INLINE_OPT void VattentionSystolicArray___024root___ico_sequent__TOP__0(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___ico_sequent__TOP__0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelfRef.attentionSystolicArray__DOT__doProcess_d 
        = (((~ (IData)(vlSelfRef.attentionSystolicArray__DOT__weight_load_enable_q)) 
            & (IData)(vlSelfRef.valid_input)) | ((0x32U 
                                                  != vlSelfRef.attentionSystolicArray__DOT__counter_q) 
                                                 & (IData)(vlSelfRef.attentionSystolicArray__DOT__doProcess_q)));
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

#ifdef VL_DEBUG
VL_ATTR_COLD void VattentionSystolicArray___024root___dump_triggers__act(VattentionSystolicArray___024root* vlSelf);
#endif  // VL_DEBUG

void VattentionSystolicArray___024root___eval_triggers__act(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_triggers__act\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    vlSelfRef.__VactTriggered.set(0U, ((IData)(vlSelfRef.clk) 
                                       & (~ (IData)(vlSelfRef.__Vtrigprevexpr___TOP__clk__0))));
    vlSelfRef.__Vtrigprevexpr___TOP__clk__0 = vlSelfRef.clk;
#ifdef VL_DEBUG
    if (VL_UNLIKELY(vlSymsp->_vm_contextp__->debug())) {
        VattentionSystolicArray___024root___dump_triggers__act(vlSelf);
    }
#endif
}

void VattentionSystolicArray___024root___nba_sequent__TOP__0(VattentionSystolicArray___024root* vlSelf);
void VattentionSystolicArray_systolic_array__N4_M4___nba_sequent__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0(VattentionSystolicArray_systolic_array__N4_M4* vlSelf);
void VattentionSystolicArray_systolic_array__N4_M4___nba_sequent__TOP__attentionSystolicArray__DOT__V_mult_systolicArray__0(VattentionSystolicArray_systolic_array__N4_M4* vlSelf);
void VattentionSystolicArray___024root___nba_sequent__TOP__1(VattentionSystolicArray___024root* vlSelf);

void VattentionSystolicArray___024root___eval_nba(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___eval_nba\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Body
    if ((1ULL & vlSelfRef.__VnbaTriggered.word(0U))) {
        VattentionSystolicArray___024root___nba_sequent__TOP__0(vlSelf);
        vlSelfRef.__Vm_traceActivity[2U] = 1U;
        VattentionSystolicArray_systolic_array__N4_M4___nba_sequent__TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray__0((&vlSymsp->TOP__attentionSystolicArray__DOT__K_mult_Q_systolicArray));
        VattentionSystolicArray_systolic_array__N4_M4___nba_sequent__TOP__attentionSystolicArray__DOT__V_mult_systolicArray__0((&vlSymsp->TOP__attentionSystolicArray__DOT__V_mult_systolicArray));
        VattentionSystolicArray___024root___nba_sequent__TOP__1(vlSelf);
    }
}

VL_INLINE_OPT void VattentionSystolicArray___024root___nba_sequent__TOP__0(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___nba_sequent__TOP__0\n"); );
    auto &vlSelfRef = std::ref(*vlSelf).get();
    // Init
    CData/*1:0*/ __Vdly__attentionSystolicArray__DOT__weight_load_counter;
    __Vdly__attentionSystolicArray__DOT__weight_load_counter = 0;
    CData/*0:0*/ __VdlySet__attentionSystolicArray__DOT__data_in_array__v0;
    __VdlySet__attentionSystolicArray__DOT__data_in_array__v0 = 0;
    double __VdlyVal__attentionSystolicArray__DOT__data_in_array__v4;
    __VdlyVal__attentionSystolicArray__DOT__data_in_array__v4 = 0;
    CData/*0:0*/ __VdlySet__attentionSystolicArray__DOT__data_in_array__v4;
    __VdlySet__attentionSystolicArray__DOT__data_in_array__v4 = 0;
    double __VdlyVal__attentionSystolicArray__DOT__data_in_array__v5;
    __VdlyVal__attentionSystolicArray__DOT__data_in_array__v5 = 0;
    double __VdlyVal__attentionSystolicArray__DOT__data_in_array__v6;
    __VdlyVal__attentionSystolicArray__DOT__data_in_array__v6 = 0;
    double __VdlyVal__attentionSystolicArray__DOT__data_in_array__v7;
    __VdlyVal__attentionSystolicArray__DOT__data_in_array__v7 = 0;
    double __VdlyVal__Q_mult_K__v0;
    __VdlyVal__Q_mult_K__v0 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v0;
    __VdlySet__Q_mult_K__v0 = 0;
    double __VdlyVal__Q_mult_K__v1;
    __VdlyVal__Q_mult_K__v1 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v1;
    __VdlySet__Q_mult_K__v1 = 0;
    double __VdlyVal__Q_mult_K__v2;
    __VdlyVal__Q_mult_K__v2 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v2;
    __VdlySet__Q_mult_K__v2 = 0;
    double __VdlyVal__Q_mult_K__v3;
    __VdlyVal__Q_mult_K__v3 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v3;
    __VdlySet__Q_mult_K__v3 = 0;
    double __VdlyVal__Q_mult_K__v4;
    __VdlyVal__Q_mult_K__v4 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v4;
    __VdlySet__Q_mult_K__v4 = 0;
    double __VdlyVal__Q_mult_K__v5;
    __VdlyVal__Q_mult_K__v5 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v5;
    __VdlySet__Q_mult_K__v5 = 0;
    double __VdlyVal__Q_mult_K__v6;
    __VdlyVal__Q_mult_K__v6 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v6;
    __VdlySet__Q_mult_K__v6 = 0;
    double __VdlyVal__Q_mult_K__v7;
    __VdlyVal__Q_mult_K__v7 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v7;
    __VdlySet__Q_mult_K__v7 = 0;
    double __VdlyVal__Q_mult_K__v8;
    __VdlyVal__Q_mult_K__v8 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v8;
    __VdlySet__Q_mult_K__v8 = 0;
    double __VdlyVal__Q_mult_K__v9;
    __VdlyVal__Q_mult_K__v9 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v9;
    __VdlySet__Q_mult_K__v9 = 0;
    double __VdlyVal__Q_mult_K__v10;
    __VdlyVal__Q_mult_K__v10 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v10;
    __VdlySet__Q_mult_K__v10 = 0;
    double __VdlyVal__Q_mult_K__v11;
    __VdlyVal__Q_mult_K__v11 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v11;
    __VdlySet__Q_mult_K__v11 = 0;
    double __VdlyVal__Q_mult_K__v12;
    __VdlyVal__Q_mult_K__v12 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v12;
    __VdlySet__Q_mult_K__v12 = 0;
    double __VdlyVal__Q_mult_K__v13;
    __VdlyVal__Q_mult_K__v13 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v13;
    __VdlySet__Q_mult_K__v13 = 0;
    double __VdlyVal__Q_mult_K__v14;
    __VdlyVal__Q_mult_K__v14 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v14;
    __VdlySet__Q_mult_K__v14 = 0;
    double __VdlyVal__Q_mult_K__v15;
    __VdlyVal__Q_mult_K__v15 = 0;
    CData/*0:0*/ __VdlySet__Q_mult_K__v15;
    __VdlySet__Q_mult_K__v15 = 0;
    double __VdlyVal__exp_Q_mult_K__v0;
    __VdlyVal__exp_Q_mult_K__v0 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v0;
    __VdlySet__exp_Q_mult_K__v0 = 0;
    double __VdlyVal__exp_Q_mult_K__v1;
    __VdlyVal__exp_Q_mult_K__v1 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v1;
    __VdlySet__exp_Q_mult_K__v1 = 0;
    double __VdlyVal__exp_Q_mult_K__v2;
    __VdlyVal__exp_Q_mult_K__v2 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v2;
    __VdlySet__exp_Q_mult_K__v2 = 0;
    double __VdlyVal__exp_Q_mult_K__v3;
    __VdlyVal__exp_Q_mult_K__v3 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v3;
    __VdlySet__exp_Q_mult_K__v3 = 0;
    double __VdlyVal__exp_Q_mult_K__v4;
    __VdlyVal__exp_Q_mult_K__v4 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v4;
    __VdlySet__exp_Q_mult_K__v4 = 0;
    double __VdlyVal__exp_Q_mult_K__v5;
    __VdlyVal__exp_Q_mult_K__v5 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v5;
    __VdlySet__exp_Q_mult_K__v5 = 0;
    double __VdlyVal__exp_Q_mult_K__v6;
    __VdlyVal__exp_Q_mult_K__v6 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v6;
    __VdlySet__exp_Q_mult_K__v6 = 0;
    double __VdlyVal__exp_Q_mult_K__v7;
    __VdlyVal__exp_Q_mult_K__v7 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v7;
    __VdlySet__exp_Q_mult_K__v7 = 0;
    double __VdlyVal__exp_Q_mult_K__v8;
    __VdlyVal__exp_Q_mult_K__v8 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v8;
    __VdlySet__exp_Q_mult_K__v8 = 0;
    double __VdlyVal__exp_Q_mult_K__v9;
    __VdlyVal__exp_Q_mult_K__v9 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v9;
    __VdlySet__exp_Q_mult_K__v9 = 0;
    double __VdlyVal__exp_Q_mult_K__v10;
    __VdlyVal__exp_Q_mult_K__v10 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v10;
    __VdlySet__exp_Q_mult_K__v10 = 0;
    double __VdlyVal__exp_Q_mult_K__v11;
    __VdlyVal__exp_Q_mult_K__v11 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v11;
    __VdlySet__exp_Q_mult_K__v11 = 0;
    double __VdlyVal__exp_Q_mult_K__v12;
    __VdlyVal__exp_Q_mult_K__v12 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v12;
    __VdlySet__exp_Q_mult_K__v12 = 0;
    double __VdlyVal__exp_Q_mult_K__v13;
    __VdlyVal__exp_Q_mult_K__v13 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v13;
    __VdlySet__exp_Q_mult_K__v13 = 0;
    double __VdlyVal__exp_Q_mult_K__v14;
    __VdlyVal__exp_Q_mult_K__v14 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v14;
    __VdlySet__exp_Q_mult_K__v14 = 0;
    double __VdlyVal__exp_Q_mult_K__v15;
    __VdlyVal__exp_Q_mult_K__v15 = 0;
    CData/*0:0*/ __VdlySet__exp_Q_mult_K__v15;
    __VdlySet__exp_Q_mult_K__v15 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v0;
    __VdlyVal__exp_K_mult_Q_mult_V__v0 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v0;
    __VdlySet__exp_K_mult_Q_mult_V__v0 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v1;
    __VdlyVal__exp_K_mult_Q_mult_V__v1 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v1;
    __VdlySet__exp_K_mult_Q_mult_V__v1 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v2;
    __VdlyVal__exp_K_mult_Q_mult_V__v2 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v2;
    __VdlySet__exp_K_mult_Q_mult_V__v2 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v3;
    __VdlyVal__exp_K_mult_Q_mult_V__v3 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v3;
    __VdlySet__exp_K_mult_Q_mult_V__v3 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v4;
    __VdlyVal__exp_K_mult_Q_mult_V__v4 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v4;
    __VdlySet__exp_K_mult_Q_mult_V__v4 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v5;
    __VdlyVal__exp_K_mult_Q_mult_V__v5 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v5;
    __VdlySet__exp_K_mult_Q_mult_V__v5 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v6;
    __VdlyVal__exp_K_mult_Q_mult_V__v6 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v6;
    __VdlySet__exp_K_mult_Q_mult_V__v6 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v7;
    __VdlyVal__exp_K_mult_Q_mult_V__v7 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v7;
    __VdlySet__exp_K_mult_Q_mult_V__v7 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v8;
    __VdlyVal__exp_K_mult_Q_mult_V__v8 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v8;
    __VdlySet__exp_K_mult_Q_mult_V__v8 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v9;
    __VdlyVal__exp_K_mult_Q_mult_V__v9 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v9;
    __VdlySet__exp_K_mult_Q_mult_V__v9 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v10;
    __VdlyVal__exp_K_mult_Q_mult_V__v10 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v10;
    __VdlySet__exp_K_mult_Q_mult_V__v10 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v11;
    __VdlyVal__exp_K_mult_Q_mult_V__v11 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v11;
    __VdlySet__exp_K_mult_Q_mult_V__v11 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v12;
    __VdlyVal__exp_K_mult_Q_mult_V__v12 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v12;
    __VdlySet__exp_K_mult_Q_mult_V__v12 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v13;
    __VdlyVal__exp_K_mult_Q_mult_V__v13 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v13;
    __VdlySet__exp_K_mult_Q_mult_V__v13 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v14;
    __VdlyVal__exp_K_mult_Q_mult_V__v14 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v14;
    __VdlySet__exp_K_mult_Q_mult_V__v14 = 0;
    double __VdlyVal__exp_K_mult_Q_mult_V__v15;
    __VdlyVal__exp_K_mult_Q_mult_V__v15 = 0;
    CData/*0:0*/ __VdlySet__exp_K_mult_Q_mult_V__v15;
    __VdlySet__exp_K_mult_Q_mult_V__v15 = 0;
    double __VdlyVal__attention__v0;
    __VdlyVal__attention__v0 = 0;
    CData/*0:0*/ __VdlySet__attention__v0;
    __VdlySet__attention__v0 = 0;
    double __VdlyVal__attention__v1;
    __VdlyVal__attention__v1 = 0;
    CData/*0:0*/ __VdlySet__attention__v1;
    __VdlySet__attention__v1 = 0;
    double __VdlyVal__attention__v2;
    __VdlyVal__attention__v2 = 0;
    CData/*0:0*/ __VdlySet__attention__v2;
    __VdlySet__attention__v2 = 0;
    double __VdlyVal__attention__v3;
    __VdlyVal__attention__v3 = 0;
    CData/*0:0*/ __VdlySet__attention__v3;
    __VdlySet__attention__v3 = 0;
    double __VdlyVal__attention__v4;
    __VdlyVal__attention__v4 = 0;
    CData/*0:0*/ __VdlySet__attention__v4;
    __VdlySet__attention__v4 = 0;
    double __VdlyVal__attention__v5;
    __VdlyVal__attention__v5 = 0;
    CData/*0:0*/ __VdlySet__attention__v5;
    __VdlySet__attention__v5 = 0;
    double __VdlyVal__attention__v6;
    __VdlyVal__attention__v6 = 0;
    CData/*0:0*/ __VdlySet__attention__v6;
    __VdlySet__attention__v6 = 0;
    double __VdlyVal__attention__v7;
    __VdlyVal__attention__v7 = 0;
    CData/*0:0*/ __VdlySet__attention__v7;
    __VdlySet__attention__v7 = 0;
    double __VdlyVal__attention__v8;
    __VdlyVal__attention__v8 = 0;
    CData/*0:0*/ __VdlySet__attention__v8;
    __VdlySet__attention__v8 = 0;
    double __VdlyVal__attention__v9;
    __VdlyVal__attention__v9 = 0;
    CData/*0:0*/ __VdlySet__attention__v9;
    __VdlySet__attention__v9 = 0;
    double __VdlyVal__attention__v10;
    __VdlyVal__attention__v10 = 0;
    CData/*0:0*/ __VdlySet__attention__v10;
    __VdlySet__attention__v10 = 0;
    double __VdlyVal__attention__v11;
    __VdlyVal__attention__v11 = 0;
    CData/*0:0*/ __VdlySet__attention__v11;
    __VdlySet__attention__v11 = 0;
    double __VdlyVal__attention__v12;
    __VdlyVal__attention__v12 = 0;
    CData/*0:0*/ __VdlySet__attention__v12;
    __VdlySet__attention__v12 = 0;
    double __VdlyVal__attention__v13;
    __VdlyVal__attention__v13 = 0;
    CData/*0:0*/ __VdlySet__attention__v13;
    __VdlySet__attention__v13 = 0;
    double __VdlyVal__attention__v14;
    __VdlyVal__attention__v14 = 0;
    CData/*0:0*/ __VdlySet__attention__v14;
    __VdlySet__attention__v14 = 0;
    double __VdlyVal__attention__v15;
    __VdlyVal__attention__v15 = 0;
    CData/*0:0*/ __VdlySet__attention__v15;
    __VdlySet__attention__v15 = 0;
    // Body
    __Vdly__attentionSystolicArray__DOT__weight_load_counter 
        = vlSelfRef.attentionSystolicArray__DOT__weight_load_counter;
    vlSelfRef.__Vdly__attentionSystolicArray__DOT__weight_load_enable_q 
        = vlSelfRef.attentionSystolicArray__DOT__weight_load_enable_q;
    __VdlySet__attention__v0 = 0U;
    __VdlySet__attention__v1 = 0U;
    __VdlySet__attention__v2 = 0U;
    __VdlySet__attention__v3 = 0U;
    __VdlySet__attention__v4 = 0U;
    __VdlySet__attention__v5 = 0U;
    __VdlySet__attention__v6 = 0U;
    __VdlySet__attention__v7 = 0U;
    __VdlySet__attention__v8 = 0U;
    __VdlySet__attention__v9 = 0U;
    __VdlySet__attention__v10 = 0U;
    __VdlySet__attention__v11 = 0U;
    __VdlySet__attention__v12 = 0U;
    __VdlySet__attention__v13 = 0U;
    __VdlySet__attention__v14 = 0U;
    __VdlySet__attention__v15 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v0 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v1 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v2 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v3 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v4 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v5 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v6 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v7 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v8 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v9 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v10 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v11 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v12 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v13 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v14 = 0U;
    __VdlySet__exp_K_mult_Q_mult_V__v15 = 0U;
    __VdlySet__exp_Q_mult_K__v0 = 0U;
    __VdlySet__exp_Q_mult_K__v1 = 0U;
    __VdlySet__exp_Q_mult_K__v2 = 0U;
    __VdlySet__exp_Q_mult_K__v3 = 0U;
    __VdlySet__exp_Q_mult_K__v4 = 0U;
    __VdlySet__exp_Q_mult_K__v5 = 0U;
    __VdlySet__exp_Q_mult_K__v6 = 0U;
    __VdlySet__exp_Q_mult_K__v7 = 0U;
    __VdlySet__exp_Q_mult_K__v8 = 0U;
    __VdlySet__exp_Q_mult_K__v9 = 0U;
    __VdlySet__exp_Q_mult_K__v10 = 0U;
    __VdlySet__exp_Q_mult_K__v11 = 0U;
    __VdlySet__exp_Q_mult_K__v12 = 0U;
    __VdlySet__exp_Q_mult_K__v13 = 0U;
    __VdlySet__exp_Q_mult_K__v14 = 0U;
    __VdlySet__exp_Q_mult_K__v15 = 0U;
    __VdlySet__attentionSystolicArray__DOT__data_in_array__v0 = 0U;
    __VdlySet__attentionSystolicArray__DOT__data_in_array__v4 = 0U;
    __VdlySet__Q_mult_K__v0 = 0U;
    __VdlySet__Q_mult_K__v1 = 0U;
    __VdlySet__Q_mult_K__v2 = 0U;
    __VdlySet__Q_mult_K__v3 = 0U;
    __VdlySet__Q_mult_K__v4 = 0U;
    __VdlySet__Q_mult_K__v5 = 0U;
    __VdlySet__Q_mult_K__v6 = 0U;
    __VdlySet__Q_mult_K__v7 = 0U;
    __VdlySet__Q_mult_K__v8 = 0U;
    __VdlySet__Q_mult_K__v9 = 0U;
    __VdlySet__Q_mult_K__v10 = 0U;
    __VdlySet__Q_mult_K__v11 = 0U;
    __VdlySet__Q_mult_K__v12 = 0U;
    __VdlySet__Q_mult_K__v13 = 0U;
    __VdlySet__Q_mult_K__v14 = 0U;
    __VdlySet__Q_mult_K__v15 = 0U;
    vlSelfRef.attentionSystolicArray__DOT__unnamedblk9__DOT__unnamedblk10__DOT__j = 4U;
    vlSelfRef.attentionSystolicArray__DOT__unnamedblk9__DOT__i = 4U;
    vlSelfRef.attentionSystolicArray__DOT__unnamedblk11__DOT__unnamedblk12__DOT__j = 4U;
    vlSelfRef.attentionSystolicArray__DOT__unnamedblk11__DOT__i = 4U;
    vlSelfRef.attentionSystolicArray__DOT__unnamedblk13__DOT__unnamedblk14__DOT__j = 4U;
    vlSelfRef.attentionSystolicArray__DOT__unnamedblk13__DOT__i = 4U;
    vlSelfRef.attentionSystolicArray__DOT__unnamedblk15__DOT__unnamedblk16__DOT__j = 4U;
    vlSelfRef.attentionSystolicArray__DOT__unnamedblk15__DOT__i = 4U;
    if ((1U & (~ (IData)(vlSelfRef.reset)))) {
        if (vlSelfRef.attentionSystolicArray__DOT__doProcess_d) {
            vlSelfRef.attentionSystolicArray__DOT__unnamedblk2__DOT__j = 4U;
        }
    }
    if ((0x29U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__attention__v0 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [3U];
        __VdlySet__attention__v0 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v1 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [3U];
        __VdlySet__exp_K_mult_Q_mult_V__v1 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v4 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [2U];
        __VdlySet__exp_K_mult_Q_mult_V__v4 = 1U;
        __VdlyVal__exp_Q_mult_K__v11 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [1U];
        __VdlySet__exp_Q_mult_K__v11 = 1U;
        __VdlyVal__exp_Q_mult_K__v14 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [0U];
        __VdlySet__exp_Q_mult_K__v14 = 1U;
    }
    if ((0x2aU == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__attention__v1 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [3U];
        __VdlySet__attention__v1 = 1U;
        __VdlyVal__attention__v4 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [2U];
        __VdlySet__attention__v4 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v2 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [3U];
        __VdlySet__exp_K_mult_Q_mult_V__v2 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v5 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [2U];
        __VdlySet__exp_K_mult_Q_mult_V__v5 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v8 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [1U];
        __VdlySet__exp_K_mult_Q_mult_V__v8 = 1U;
        __VdlyVal__exp_Q_mult_K__v15 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [0U];
        __VdlySet__exp_Q_mult_K__v15 = 1U;
    }
    if ((0x2bU == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__attention__v2 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [3U];
        __VdlySet__attention__v2 = 1U;
        __VdlyVal__attention__v5 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [2U];
        __VdlySet__attention__v5 = 1U;
        __VdlyVal__attention__v8 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [1U];
        __VdlySet__attention__v8 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v3 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [3U];
        __VdlySet__exp_K_mult_Q_mult_V__v3 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v6 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [2U];
        __VdlySet__exp_K_mult_Q_mult_V__v6 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v9 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [1U];
        __VdlySet__exp_K_mult_Q_mult_V__v9 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v12 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [0U];
        __VdlySet__exp_K_mult_Q_mult_V__v12 = 1U;
    }
    if ((0x2cU == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__attention__v3 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [3U];
        __VdlySet__attention__v3 = 1U;
        __VdlyVal__attention__v6 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [2U];
        __VdlySet__attention__v6 = 1U;
        __VdlyVal__attention__v9 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [1U];
        __VdlySet__attention__v9 = 1U;
        __VdlyVal__attention__v12 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [0U];
        __VdlySet__attention__v12 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v7 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [2U];
        __VdlySet__exp_K_mult_Q_mult_V__v7 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v10 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [1U];
        __VdlySet__exp_K_mult_Q_mult_V__v10 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v13 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [0U];
        __VdlySet__exp_K_mult_Q_mult_V__v13 = 1U;
    }
    if ((0x2dU == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__attention__v7 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [2U];
        __VdlySet__attention__v7 = 1U;
        __VdlyVal__attention__v10 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [1U];
        __VdlySet__attention__v10 = 1U;
        __VdlyVal__attention__v13 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [0U];
        __VdlySet__attention__v13 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v11 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [1U];
        __VdlySet__exp_K_mult_Q_mult_V__v11 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v14 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [0U];
        __VdlySet__exp_K_mult_Q_mult_V__v14 = 1U;
    }
    if ((0x2eU == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__attention__v11 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [1U];
        __VdlySet__attention__v11 = 1U;
        __VdlyVal__attention__v14 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [0U];
        __VdlySet__attention__v14 = 1U;
        __VdlyVal__exp_K_mult_Q_mult_V__v15 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [0U];
        __VdlySet__exp_K_mult_Q_mult_V__v15 = 1U;
    }
    if ((0x2fU == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__attention__v15 = vlSelfRef.attentionSystolicArray__DOT__attention_result
            [0U];
        __VdlySet__attention__v15 = 1U;
    }
    if ((0x28U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__exp_K_mult_Q_mult_V__v0 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
            [3U];
        __VdlySet__exp_K_mult_Q_mult_V__v0 = 1U;
        __VdlyVal__exp_Q_mult_K__v7 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [2U];
        __VdlySet__exp_Q_mult_K__v7 = 1U;
        __VdlyVal__exp_Q_mult_K__v10 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [1U];
        __VdlySet__exp_Q_mult_K__v10 = 1U;
        __VdlyVal__exp_Q_mult_K__v13 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [0U];
        __VdlySet__exp_Q_mult_K__v13 = 1U;
    }
    if ((0x24U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__exp_Q_mult_K__v0 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [3U];
        __VdlySet__exp_Q_mult_K__v0 = 1U;
    }
    if ((0x25U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__exp_Q_mult_K__v1 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [3U];
        __VdlySet__exp_Q_mult_K__v1 = 1U;
        __VdlyVal__exp_Q_mult_K__v4 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [2U];
        __VdlySet__exp_Q_mult_K__v4 = 1U;
    }
    if ((0x26U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__exp_Q_mult_K__v2 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [3U];
        __VdlySet__exp_Q_mult_K__v2 = 1U;
        __VdlyVal__exp_Q_mult_K__v5 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [2U];
        __VdlySet__exp_Q_mult_K__v5 = 1U;
        __VdlyVal__exp_Q_mult_K__v8 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [1U];
        __VdlySet__exp_Q_mult_K__v8 = 1U;
    }
    if ((0x27U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__exp_Q_mult_K__v3 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [3U];
        __VdlySet__exp_Q_mult_K__v3 = 1U;
        __VdlyVal__exp_Q_mult_K__v6 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [2U];
        __VdlySet__exp_Q_mult_K__v6 = 1U;
        __VdlyVal__exp_Q_mult_K__v9 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [1U];
        __VdlySet__exp_Q_mult_K__v9 = 1U;
        __VdlyVal__exp_Q_mult_K__v12 = vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_result
            [0U];
        __VdlySet__exp_Q_mult_K__v12 = 1U;
    }
    if ((4U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__Q_mult_K__v0 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [3U];
        __VdlySet__Q_mult_K__v0 = 1U;
    }
    if ((5U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__Q_mult_K__v1 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [3U];
        __VdlySet__Q_mult_K__v1 = 1U;
        __VdlyVal__Q_mult_K__v4 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [2U];
        __VdlySet__Q_mult_K__v4 = 1U;
    }
    if ((6U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__Q_mult_K__v2 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [3U];
        __VdlySet__Q_mult_K__v2 = 1U;
        __VdlyVal__Q_mult_K__v5 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [2U];
        __VdlySet__Q_mult_K__v5 = 1U;
        __VdlyVal__Q_mult_K__v8 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [1U];
        __VdlySet__Q_mult_K__v8 = 1U;
    }
    if ((7U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__Q_mult_K__v3 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [3U];
        __VdlySet__Q_mult_K__v3 = 1U;
        __VdlyVal__Q_mult_K__v6 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [2U];
        __VdlySet__Q_mult_K__v6 = 1U;
        __VdlyVal__Q_mult_K__v9 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [1U];
        __VdlySet__Q_mult_K__v9 = 1U;
        __VdlyVal__Q_mult_K__v12 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [0U];
        __VdlySet__Q_mult_K__v12 = 1U;
    }
    if ((8U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__Q_mult_K__v7 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [2U];
        __VdlySet__Q_mult_K__v7 = 1U;
        __VdlyVal__Q_mult_K__v10 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [1U];
        __VdlySet__Q_mult_K__v10 = 1U;
        __VdlyVal__Q_mult_K__v13 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [0U];
        __VdlySet__Q_mult_K__v13 = 1U;
    }
    if ((9U == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__Q_mult_K__v11 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [1U];
        __VdlySet__Q_mult_K__v11 = 1U;
        __VdlyVal__Q_mult_K__v14 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [0U];
        __VdlySet__Q_mult_K__v14 = 1U;
    }
    if ((0xaU == vlSelfRef.attentionSystolicArray__DOT__counter_q)) {
        __VdlyVal__Q_mult_K__v15 = vlSelfRef.attentionSystolicArray__DOT__Q_mult_K_result
            [0U];
        __VdlySet__Q_mult_K__v15 = 1U;
    }
    vlSelfRef.attentionSystolicArray__DOT__valid_result_q 
        = ((~ (IData)(vlSelfRef.reset)) & (0x31U == vlSelfRef.attentionSystolicArray__DOT__counter_q));
    if (vlSelfRef.reset) {
        vlSelfRef.attentionSystolicArray__DOT__unnamedblk1__DOT__j = 4U;
        vlSelfRef.__Vdly__attentionSystolicArray__DOT__weight_load_enable_q = 1U;
        __Vdly__attentionSystolicArray__DOT__weight_load_counter = 0U;
        __VdlySet__attentionSystolicArray__DOT__data_in_array__v0 = 1U;
        vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__dummy_pe__data_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__norm_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__norm_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__norm_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__norm_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__x_normed_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__x_normed_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__x_normed_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__x_normed_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__x_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_sum_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_sum_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_sum_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_sum_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_out = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__counter_q = 0U;
    } else {
        if (vlSelfRef.attentionSystolicArray__DOT__weight_load_enable_q) {
            __Vdly__attentionSystolicArray__DOT__weight_load_counter 
                = (3U & ((IData)(1U) + (IData)(vlSelfRef.attentionSystolicArray__DOT__weight_load_counter)));
            if ((2U == (IData)(vlSelfRef.attentionSystolicArray__DOT__weight_load_counter))) {
                vlSelfRef.__Vdly__attentionSystolicArray__DOT__weight_load_enable_q = 0U;
            }
        }
        if (vlSelfRef.attentionSystolicArray__DOT__doProcess_d) {
            __VdlyVal__attentionSystolicArray__DOT__data_in_array__v4 
                = ((VL_LTES_III(32, 0U, vlSelfRef.attentionSystolicArray__DOT__counter_d) 
                    & (4U > vlSelfRef.attentionSystolicArray__DOT__counter_d))
                    ? vlSelfRef.Q_matrix[(3U & vlSelfRef.attentionSystolicArray__DOT__counter_d)]
                   [0U] : 0.0);
            __VdlySet__attentionSystolicArray__DOT__data_in_array__v4 = 1U;
            __VdlyVal__attentionSystolicArray__DOT__data_in_array__v5 
                = ((VL_LTES_III(32, 0U, (vlSelfRef.attentionSystolicArray__DOT__counter_d 
                                         - (IData)(1U))) 
                    & (4U > (vlSelfRef.attentionSystolicArray__DOT__counter_d 
                             - (IData)(1U)))) ? vlSelfRef.Q_matrix
                   [(3U & (vlSelfRef.attentionSystolicArray__DOT__counter_d 
                           - (IData)(1U)))][1U] : 0.0);
            __VdlyVal__attentionSystolicArray__DOT__data_in_array__v6 
                = ((VL_LTES_III(32, 0U, (vlSelfRef.attentionSystolicArray__DOT__counter_d 
                                         - (IData)(2U))) 
                    & (4U > (vlSelfRef.attentionSystolicArray__DOT__counter_d 
                             - (IData)(2U)))) ? vlSelfRef.Q_matrix
                   [(3U & (vlSelfRef.attentionSystolicArray__DOT__counter_d 
                           - (IData)(2U)))][2U] : 0.0);
            __VdlyVal__attentionSystolicArray__DOT__data_in_array__v7 
                = ((VL_LTES_III(32, 0U, (vlSelfRef.attentionSystolicArray__DOT__counter_d 
                                         - (IData)(3U))) 
                    & (4U > (vlSelfRef.attentionSystolicArray__DOT__counter_d 
                             - (IData)(3U)))) ? vlSelfRef.Q_matrix
                   [(3U & (vlSelfRef.attentionSystolicArray__DOT__counter_d 
                           - (IData)(3U)))][3U] : 0.0);
        }
        if (vlSelfRef.attentionSystolicArray__DOT__doProcess_q) {
            vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__dummy_pe__data_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                [0U];
            vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__norm_out 
                = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                [1U];
            vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__norm_out 
                = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                [2U];
            vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__norm_out 
                = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                [3U];
            vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__norm_out 
                = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                [4U];
            vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__x_normed_out 
                = (vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                   [0U] / vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                   [1U]);
            vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__x_normed_out 
                = (vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                   [1U] / vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                   [2U]);
            vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__x_normed_out 
                = (vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                   [2U] / vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                   [3U]);
            vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__x_normed_out 
                = (vlSelfRef.attentionSystolicArray__DOT__exp_Q_mult_K_mult_V_result
                   [3U] / vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT__norm_wire
                   [4U]);
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__x_pow_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__exp_term;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__taylor_approx_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__DOT__taylor_sum;
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__0__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__1__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][1U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__2__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][2U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__3__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][3U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__4__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][4U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__5__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][5U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__6__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][6U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__7__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][7U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__8__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][8U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__9__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][9U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__10__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0xaU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__11__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0xbU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__12__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0xcU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__13__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0xdU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__14__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0xeU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__15__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0xfU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__16__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x10U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__17__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x11U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__18__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x12U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__19__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x13U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__20__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x14U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__21__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x15U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__22__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x16U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__23__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x17U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__24__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x18U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__25__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x19U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__26__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x1aU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__27__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x1bU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__28__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x1cU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__29__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x1dU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__col__BRA__30__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [0U][0x1eU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__0__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__1__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][1U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__2__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][2U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__3__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][3U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__4__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][4U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__5__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][5U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__6__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][6U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__7__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][7U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__8__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][8U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__9__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][9U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__10__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0xaU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__11__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0xbU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__12__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0xcU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__13__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0xdU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__14__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0xeU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__15__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0xfU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__16__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x10U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__17__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x11U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__18__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x12U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__19__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x13U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__20__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x14U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__21__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x15U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__22__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x16U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__23__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x17U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__24__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x18U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__25__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x19U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__26__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x1aU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__27__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x1bU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__28__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x1cU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__29__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x1dU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__col__BRA__30__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [1U][0x1eU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__0__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__1__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][1U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__2__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][2U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__3__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][3U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__4__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][4U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__5__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][5U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__6__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][6U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__7__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][7U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__8__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][8U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__9__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][9U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__10__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0xaU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__11__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0xbU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__12__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0xcU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__13__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0xdU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__14__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0xeU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__15__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0xfU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__16__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x10U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__17__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x11U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__18__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x12U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__19__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x13U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__20__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x14U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__21__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x15U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__22__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x16U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__23__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x17U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__24__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x18U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__25__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x19U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__26__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x1aU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__27__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x1bU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__28__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x1cU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__29__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x1dU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__col__BRA__30__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [2U][0x1eU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__0__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__1__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][1U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__2__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][2U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__3__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][3U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__4__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][4U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__5__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][5U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__6__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][6U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__7__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][7U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__8__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][8U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__9__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][9U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__10__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0xaU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__11__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0xbU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__12__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0xcU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__13__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0xdU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__14__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0xeU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__15__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0xfU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__16__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x10U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__17__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x11U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__18__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x12U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__19__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x13U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__20__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x14U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__21__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x15U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__22__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x16U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__23__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x17U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__24__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x18U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__25__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x19U];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__26__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x1aU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__27__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x1bU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__28__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x1cU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__29__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x1dU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__col__BRA__30__KET____DOT__pe__x_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__x_wire
                [3U][0x1eU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_sum_out 
                = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                   [0U][0x1fU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                   [1U]);
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_sum_out 
                = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                   [1U][0x1fU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                   [2U]);
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_sum_out 
                = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                   [2U][0x1fU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                   [3U]);
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_sum_out 
                = (vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                   [3U][0x1fU] + vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire
                   [4U]);
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                [0U][0x1fU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                [1U][0x1fU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                [2U][0x1fU];
            vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_out 
                = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire
                [3U][0x1fU];
        }
        vlSelfRef.attentionSystolicArray__DOT__counter_q 
            = vlSelfRef.attentionSystolicArray__DOT__counter_d;
    }
    vlSelfRef.attentionSystolicArray__DOT__weight_load_counter 
        = __Vdly__attentionSystolicArray__DOT__weight_load_counter;
    if (__VdlySet__attention__v0) {
        vlSelfRef.attention[0U][0U] = __VdlyVal__attention__v0;
    }
    if (__VdlySet__attention__v1) {
        vlSelfRef.attention[1U][0U] = __VdlyVal__attention__v1;
    }
    if (__VdlySet__attention__v2) {
        vlSelfRef.attention[2U][0U] = __VdlyVal__attention__v2;
    }
    if (__VdlySet__attention__v3) {
        vlSelfRef.attention[3U][0U] = __VdlyVal__attention__v3;
    }
    if (__VdlySet__attention__v4) {
        vlSelfRef.attention[0U][1U] = __VdlyVal__attention__v4;
    }
    if (__VdlySet__attention__v5) {
        vlSelfRef.attention[1U][1U] = __VdlyVal__attention__v5;
    }
    if (__VdlySet__attention__v6) {
        vlSelfRef.attention[2U][1U] = __VdlyVal__attention__v6;
    }
    if (__VdlySet__attention__v7) {
        vlSelfRef.attention[3U][1U] = __VdlyVal__attention__v7;
    }
    if (__VdlySet__attention__v8) {
        vlSelfRef.attention[0U][2U] = __VdlyVal__attention__v8;
    }
    if (__VdlySet__attention__v9) {
        vlSelfRef.attention[1U][2U] = __VdlyVal__attention__v9;
    }
    if (__VdlySet__attention__v10) {
        vlSelfRef.attention[2U][2U] = __VdlyVal__attention__v10;
    }
    if (__VdlySet__attention__v11) {
        vlSelfRef.attention[3U][2U] = __VdlyVal__attention__v11;
    }
    if (__VdlySet__attention__v12) {
        vlSelfRef.attention[0U][3U] = __VdlyVal__attention__v12;
    }
    if (__VdlySet__attention__v13) {
        vlSelfRef.attention[1U][3U] = __VdlyVal__attention__v13;
    }
    if (__VdlySet__attention__v14) {
        vlSelfRef.attention[2U][3U] = __VdlyVal__attention__v14;
    }
    if (__VdlySet__attention__v15) {
        vlSelfRef.attention[3U][3U] = __VdlyVal__attention__v15;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v0) {
        vlSelfRef.exp_K_mult_Q_mult_V[0U][0U] = __VdlyVal__exp_K_mult_Q_mult_V__v0;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v1) {
        vlSelfRef.exp_K_mult_Q_mult_V[1U][0U] = __VdlyVal__exp_K_mult_Q_mult_V__v1;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v2) {
        vlSelfRef.exp_K_mult_Q_mult_V[2U][0U] = __VdlyVal__exp_K_mult_Q_mult_V__v2;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v3) {
        vlSelfRef.exp_K_mult_Q_mult_V[3U][0U] = __VdlyVal__exp_K_mult_Q_mult_V__v3;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v4) {
        vlSelfRef.exp_K_mult_Q_mult_V[0U][1U] = __VdlyVal__exp_K_mult_Q_mult_V__v4;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v5) {
        vlSelfRef.exp_K_mult_Q_mult_V[1U][1U] = __VdlyVal__exp_K_mult_Q_mult_V__v5;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v6) {
        vlSelfRef.exp_K_mult_Q_mult_V[2U][1U] = __VdlyVal__exp_K_mult_Q_mult_V__v6;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v7) {
        vlSelfRef.exp_K_mult_Q_mult_V[3U][1U] = __VdlyVal__exp_K_mult_Q_mult_V__v7;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v8) {
        vlSelfRef.exp_K_mult_Q_mult_V[0U][2U] = __VdlyVal__exp_K_mult_Q_mult_V__v8;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v9) {
        vlSelfRef.exp_K_mult_Q_mult_V[1U][2U] = __VdlyVal__exp_K_mult_Q_mult_V__v9;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v10) {
        vlSelfRef.exp_K_mult_Q_mult_V[2U][2U] = __VdlyVal__exp_K_mult_Q_mult_V__v10;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v11) {
        vlSelfRef.exp_K_mult_Q_mult_V[3U][2U] = __VdlyVal__exp_K_mult_Q_mult_V__v11;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v12) {
        vlSelfRef.exp_K_mult_Q_mult_V[0U][3U] = __VdlyVal__exp_K_mult_Q_mult_V__v12;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v13) {
        vlSelfRef.exp_K_mult_Q_mult_V[1U][3U] = __VdlyVal__exp_K_mult_Q_mult_V__v13;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v14) {
        vlSelfRef.exp_K_mult_Q_mult_V[2U][3U] = __VdlyVal__exp_K_mult_Q_mult_V__v14;
    }
    if (__VdlySet__exp_K_mult_Q_mult_V__v15) {
        vlSelfRef.exp_K_mult_Q_mult_V[3U][3U] = __VdlyVal__exp_K_mult_Q_mult_V__v15;
    }
    if (__VdlySet__exp_Q_mult_K__v0) {
        vlSelfRef.exp_Q_mult_K[0U][0U] = __VdlyVal__exp_Q_mult_K__v0;
    }
    if (__VdlySet__exp_Q_mult_K__v1) {
        vlSelfRef.exp_Q_mult_K[1U][0U] = __VdlyVal__exp_Q_mult_K__v1;
    }
    if (__VdlySet__exp_Q_mult_K__v2) {
        vlSelfRef.exp_Q_mult_K[2U][0U] = __VdlyVal__exp_Q_mult_K__v2;
    }
    if (__VdlySet__exp_Q_mult_K__v3) {
        vlSelfRef.exp_Q_mult_K[3U][0U] = __VdlyVal__exp_Q_mult_K__v3;
    }
    if (__VdlySet__exp_Q_mult_K__v4) {
        vlSelfRef.exp_Q_mult_K[0U][1U] = __VdlyVal__exp_Q_mult_K__v4;
    }
    if (__VdlySet__exp_Q_mult_K__v5) {
        vlSelfRef.exp_Q_mult_K[1U][1U] = __VdlyVal__exp_Q_mult_K__v5;
    }
    if (__VdlySet__exp_Q_mult_K__v6) {
        vlSelfRef.exp_Q_mult_K[2U][1U] = __VdlyVal__exp_Q_mult_K__v6;
    }
    if (__VdlySet__exp_Q_mult_K__v7) {
        vlSelfRef.exp_Q_mult_K[3U][1U] = __VdlyVal__exp_Q_mult_K__v7;
    }
    if (__VdlySet__exp_Q_mult_K__v8) {
        vlSelfRef.exp_Q_mult_K[0U][2U] = __VdlyVal__exp_Q_mult_K__v8;
    }
    if (__VdlySet__exp_Q_mult_K__v9) {
        vlSelfRef.exp_Q_mult_K[1U][2U] = __VdlyVal__exp_Q_mult_K__v9;
    }
    if (__VdlySet__exp_Q_mult_K__v10) {
        vlSelfRef.exp_Q_mult_K[2U][2U] = __VdlyVal__exp_Q_mult_K__v10;
    }
    if (__VdlySet__exp_Q_mult_K__v11) {
        vlSelfRef.exp_Q_mult_K[3U][2U] = __VdlyVal__exp_Q_mult_K__v11;
    }
    if (__VdlySet__exp_Q_mult_K__v12) {
        vlSelfRef.exp_Q_mult_K[0U][3U] = __VdlyVal__exp_Q_mult_K__v12;
    }
    if (__VdlySet__exp_Q_mult_K__v13) {
        vlSelfRef.exp_Q_mult_K[1U][3U] = __VdlyVal__exp_Q_mult_K__v13;
    }
    if (__VdlySet__exp_Q_mult_K__v14) {
        vlSelfRef.exp_Q_mult_K[2U][3U] = __VdlyVal__exp_Q_mult_K__v14;
    }
    if (__VdlySet__exp_Q_mult_K__v15) {
        vlSelfRef.exp_Q_mult_K[3U][3U] = __VdlyVal__exp_Q_mult_K__v15;
    }
    if (__VdlySet__attentionSystolicArray__DOT__data_in_array__v0) {
        vlSelfRef.attentionSystolicArray__DOT__data_in_array[0U] = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__data_in_array[1U] = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__data_in_array[2U] = 0.0;
        vlSelfRef.attentionSystolicArray__DOT__data_in_array[3U] = 0.0;
    }
    if (__VdlySet__attentionSystolicArray__DOT__data_in_array__v4) {
        vlSelfRef.attentionSystolicArray__DOT__data_in_array[0U] 
            = __VdlyVal__attentionSystolicArray__DOT__data_in_array__v4;
        vlSelfRef.attentionSystolicArray__DOT__data_in_array[1U] 
            = __VdlyVal__attentionSystolicArray__DOT__data_in_array__v5;
        vlSelfRef.attentionSystolicArray__DOT__data_in_array[2U] 
            = __VdlyVal__attentionSystolicArray__DOT__data_in_array__v6;
        vlSelfRef.attentionSystolicArray__DOT__data_in_array[3U] 
            = __VdlyVal__attentionSystolicArray__DOT__data_in_array__v7;
    }
    if (__VdlySet__Q_mult_K__v0) {
        vlSelfRef.Q_mult_K[0U][0U] = __VdlyVal__Q_mult_K__v0;
    }
    if (__VdlySet__Q_mult_K__v1) {
        vlSelfRef.Q_mult_K[1U][0U] = __VdlyVal__Q_mult_K__v1;
    }
    if (__VdlySet__Q_mult_K__v2) {
        vlSelfRef.Q_mult_K[2U][0U] = __VdlyVal__Q_mult_K__v2;
    }
    if (__VdlySet__Q_mult_K__v3) {
        vlSelfRef.Q_mult_K[3U][0U] = __VdlyVal__Q_mult_K__v3;
    }
    if (__VdlySet__Q_mult_K__v4) {
        vlSelfRef.Q_mult_K[0U][1U] = __VdlyVal__Q_mult_K__v4;
    }
    if (__VdlySet__Q_mult_K__v5) {
        vlSelfRef.Q_mult_K[1U][1U] = __VdlyVal__Q_mult_K__v5;
    }
    if (__VdlySet__Q_mult_K__v6) {
        vlSelfRef.Q_mult_K[2U][1U] = __VdlyVal__Q_mult_K__v6;
    }
    if (__VdlySet__Q_mult_K__v7) {
        vlSelfRef.Q_mult_K[3U][1U] = __VdlyVal__Q_mult_K__v7;
    }
    if (__VdlySet__Q_mult_K__v8) {
        vlSelfRef.Q_mult_K[0U][2U] = __VdlyVal__Q_mult_K__v8;
    }
    if (__VdlySet__Q_mult_K__v9) {
        vlSelfRef.Q_mult_K[1U][2U] = __VdlyVal__Q_mult_K__v9;
    }
    if (__VdlySet__Q_mult_K__v10) {
        vlSelfRef.Q_mult_K[2U][2U] = __VdlyVal__Q_mult_K__v10;
    }
    if (__VdlySet__Q_mult_K__v11) {
        vlSelfRef.Q_mult_K[3U][2U] = __VdlyVal__Q_mult_K__v11;
    }
    if (__VdlySet__Q_mult_K__v12) {
        vlSelfRef.Q_mult_K[0U][3U] = __VdlyVal__Q_mult_K__v12;
    }
    if (__VdlySet__Q_mult_K__v13) {
        vlSelfRef.Q_mult_K[1U][3U] = __VdlyVal__Q_mult_K__v13;
    }
    if (__VdlySet__Q_mult_K__v14) {
        vlSelfRef.Q_mult_K[2U][3U] = __VdlyVal__Q_mult_K__v14;
    }
    if (__VdlySet__Q_mult_K__v15) {
        vlSelfRef.Q_mult_K[3U][3U] = __VdlyVal__Q_mult_K__v15;
    }
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
    vlSelfRef.valid_result = vlSelfRef.attentionSystolicArray__DOT__valid_result_q;
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
    vlSelfRef.attentionSystolicArray__DOT__attention_result[0U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe__x_normed_out;
    vlSelfRef.attentionSystolicArray__DOT__attention_result[1U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe__x_normed_out;
    vlSelfRef.attentionSystolicArray__DOT__attention_result[2U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe__x_normed_out;
    vlSelfRef.attentionSystolicArray__DOT__attention_result[3U] 
        = vlSelfRef.attentionSystolicArray__DOT__softmax_norm_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe__x_normed_out;
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
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[0U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_sum_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[1U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_sum_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[2U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_sum_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__accumulated_sum_wire[3U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_sum_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[0U][0x20U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__0__KET____DOT__pe_acc__exp_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[1U][0x20U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__1__KET____DOT__pe_acc__exp_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[2U][0x20U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__2__KET____DOT__pe_acc__exp_out;
    vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT__taylor_approx_wire[3U][0x20U] 
        = vlSelfRef.attentionSystolicArray__DOT__exp_systolic_array__DOT____Vcellout__row__BRA__3__KET____DOT__pe_acc__exp_out;
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

VL_INLINE_OPT void VattentionSystolicArray___024root___nba_sequent__TOP__1(VattentionSystolicArray___024root* vlSelf) {
    (void)vlSelf;  // Prevent unused variable warning
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    VL_DEBUG_IF(VL_DBG_MSGF("+    VattentionSystolicArray___024root___nba_sequent__TOP__1\n"); );
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
    vlSelfRef.attentionSystolicArray__DOT__weight_load_enable_q 
        = vlSelfRef.__Vdly__attentionSystolicArray__DOT__weight_load_enable_q;
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
    vlSelfRef.attentionSystolicArray__DOT__doProcess_q 
        = ((~ (IData)(vlSelfRef.reset)) & (IData)(vlSelfRef.attentionSystolicArray__DOT__doProcess_d));
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
    vlSelfRef.attentionSystolicArray__DOT__counter_d 
        = ((IData)(vlSelfRef.attentionSystolicArray__DOT__doProcess_q)
            ? ((IData)(1U) + vlSelfRef.attentionSystolicArray__DOT__counter_q)
            : 0U);
    vlSelfRef.attentionSystolicArray__DOT__doProcess_d 
        = (((~ (IData)(vlSelfRef.attentionSystolicArray__DOT__weight_load_enable_q)) 
            & (IData)(vlSelfRef.valid_input)) | ((0x32U 
                                                  != vlSelfRef.attentionSystolicArray__DOT__counter_q) 
                                                 & (IData)(vlSelfRef.attentionSystolicArray__DOT__doProcess_q)));
}
