// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Model implementation (design independent parts)

#include "VattentionSystolicArray__pch.h"
#include "verilated_vcd_c.h"

//============================================================
// Constructors

VattentionSystolicArray::VattentionSystolicArray(VerilatedContext* _vcontextp__, const char* _vcname__)
    : VerilatedModel{*_vcontextp__}
    , vlSymsp{new VattentionSystolicArray__Syms(contextp(), _vcname__, this)}
    , clk{vlSymsp->TOP.clk}
    , reset{vlSymsp->TOP.reset}
    , valid_input{vlSymsp->TOP.valid_input}
    , valid_result{vlSymsp->TOP.valid_result}
    , softmax_temp{vlSymsp->TOP.softmax_temp}
    , K_matrix{vlSymsp->TOP.K_matrix}
    , factorial_arr{vlSymsp->TOP.factorial_arr}
    , Q_matrix{vlSymsp->TOP.Q_matrix}
    , V_matrix{vlSymsp->TOP.V_matrix}
    , Q_mult_K{vlSymsp->TOP.Q_mult_K}
    , exp_Q_mult_K{vlSymsp->TOP.exp_Q_mult_K}
    , exp_K_mult_Q_mult_V{vlSymsp->TOP.exp_K_mult_Q_mult_V}
    , attention{vlSymsp->TOP.attention}
    , __PVT__attentionSystolicArray__DOT__K_mult_Q_systolicArray{vlSymsp->TOP.__PVT__attentionSystolicArray__DOT__K_mult_Q_systolicArray}
    , __PVT__attentionSystolicArray__DOT__V_mult_systolicArray{vlSymsp->TOP.__PVT__attentionSystolicArray__DOT__V_mult_systolicArray}
    , rootp{&(vlSymsp->TOP)}
{
    // Register model with the context
    contextp()->addModel(this);
    contextp()->traceBaseModelCbAdd(
        [this](VerilatedTraceBaseC* tfp, int levels, int options) { traceBaseModel(tfp, levels, options); });
}

VattentionSystolicArray::VattentionSystolicArray(const char* _vcname__)
    : VattentionSystolicArray(Verilated::threadContextp(), _vcname__)
{
}

//============================================================
// Destructor

VattentionSystolicArray::~VattentionSystolicArray() {
    delete vlSymsp;
}

//============================================================
// Evaluation function

#ifdef VL_DEBUG
void VattentionSystolicArray___024root___eval_debug_assertions(VattentionSystolicArray___024root* vlSelf);
#endif  // VL_DEBUG
void VattentionSystolicArray___024root___eval_static(VattentionSystolicArray___024root* vlSelf);
void VattentionSystolicArray___024root___eval_initial(VattentionSystolicArray___024root* vlSelf);
void VattentionSystolicArray___024root___eval_settle(VattentionSystolicArray___024root* vlSelf);
void VattentionSystolicArray___024root___eval(VattentionSystolicArray___024root* vlSelf);

void VattentionSystolicArray::eval_step() {
    VL_DEBUG_IF(VL_DBG_MSGF("+++++TOP Evaluate VattentionSystolicArray::eval_step\n"); );
#ifdef VL_DEBUG
    // Debug assertions
    VattentionSystolicArray___024root___eval_debug_assertions(&(vlSymsp->TOP));
#endif  // VL_DEBUG
    vlSymsp->__Vm_activity = true;
    vlSymsp->__Vm_deleter.deleteAll();
    if (VL_UNLIKELY(!vlSymsp->__Vm_didInit)) {
        vlSymsp->__Vm_didInit = true;
        VL_DEBUG_IF(VL_DBG_MSGF("+ Initial\n"););
        VattentionSystolicArray___024root___eval_static(&(vlSymsp->TOP));
        VattentionSystolicArray___024root___eval_initial(&(vlSymsp->TOP));
        VattentionSystolicArray___024root___eval_settle(&(vlSymsp->TOP));
    }
    VL_DEBUG_IF(VL_DBG_MSGF("+ Eval\n"););
    VattentionSystolicArray___024root___eval(&(vlSymsp->TOP));
    // Evaluate cleanup
    Verilated::endOfEval(vlSymsp->__Vm_evalMsgQp);
}

//============================================================
// Events and timing
bool VattentionSystolicArray::eventsPending() { return false; }

uint64_t VattentionSystolicArray::nextTimeSlot() {
    VL_FATAL_MT(__FILE__, __LINE__, "", "%Error: No delays in the design");
    return 0;
}

//============================================================
// Utilities

const char* VattentionSystolicArray::name() const {
    return vlSymsp->name();
}

//============================================================
// Invoke final blocks

void VattentionSystolicArray___024root___eval_final(VattentionSystolicArray___024root* vlSelf);

VL_ATTR_COLD void VattentionSystolicArray::final() {
    VattentionSystolicArray___024root___eval_final(&(vlSymsp->TOP));
}

//============================================================
// Implementations of abstract methods from VerilatedModel

const char* VattentionSystolicArray::hierName() const { return vlSymsp->name(); }
const char* VattentionSystolicArray::modelName() const { return "VattentionSystolicArray"; }
unsigned VattentionSystolicArray::threads() const { return 1; }
void VattentionSystolicArray::prepareClone() const { contextp()->prepareClone(); }
void VattentionSystolicArray::atClone() const {
    contextp()->threadPoolpOnClone();
}
std::unique_ptr<VerilatedTraceConfig> VattentionSystolicArray::traceConfig() const {
    return std::unique_ptr<VerilatedTraceConfig>{new VerilatedTraceConfig{false, false, false}};
};

//============================================================
// Trace configuration

void VattentionSystolicArray___024root__trace_decl_types(VerilatedVcd* tracep);

void VattentionSystolicArray___024root__trace_init_top(VattentionSystolicArray___024root* vlSelf, VerilatedVcd* tracep);

VL_ATTR_COLD static void trace_init(void* voidSelf, VerilatedVcd* tracep, uint32_t code) {
    // Callback from tracep->open()
    VattentionSystolicArray___024root* const __restrict vlSelf VL_ATTR_UNUSED = static_cast<VattentionSystolicArray___024root*>(voidSelf);
    VattentionSystolicArray__Syms* const __restrict vlSymsp VL_ATTR_UNUSED = vlSelf->vlSymsp;
    if (!vlSymsp->_vm_contextp__->calcUnusedSigs()) {
        VL_FATAL_MT(__FILE__, __LINE__, __FILE__,
            "Turning on wave traces requires Verilated::traceEverOn(true) call before time 0.");
    }
    vlSymsp->__Vm_baseCode = code;
    tracep->pushPrefix(std::string{vlSymsp->name()}, VerilatedTracePrefixType::SCOPE_MODULE);
    VattentionSystolicArray___024root__trace_decl_types(tracep);
    VattentionSystolicArray___024root__trace_init_top(vlSelf, tracep);
    tracep->popPrefix();
}

VL_ATTR_COLD void VattentionSystolicArray___024root__trace_register(VattentionSystolicArray___024root* vlSelf, VerilatedVcd* tracep);

VL_ATTR_COLD void VattentionSystolicArray::traceBaseModel(VerilatedTraceBaseC* tfp, int levels, int options) {
    (void)levels; (void)options;
    VerilatedVcdC* const stfp = dynamic_cast<VerilatedVcdC*>(tfp);
    if (VL_UNLIKELY(!stfp)) {
        vl_fatal(__FILE__, __LINE__, __FILE__,"'VattentionSystolicArray::trace()' called on non-VerilatedVcdC object;"
            " use --trace-fst with VerilatedFst object, and --trace with VerilatedVcd object");
    }
    stfp->spTrace()->addModel(this);
    stfp->spTrace()->addInitCb(&trace_init, &(vlSymsp->TOP));
    VattentionSystolicArray___024root__trace_register(&(vlSymsp->TOP), stfp->spTrace());
}
