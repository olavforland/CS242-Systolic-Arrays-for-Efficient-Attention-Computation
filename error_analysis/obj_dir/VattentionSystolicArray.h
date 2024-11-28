// Verilated -*- C++ -*-
// DESCRIPTION: Verilator output: Primary model header
//
// This header should be included by all source files instantiating the design.
// The class here is then constructed to instantiate the design.
// See the Verilator manual for examples.

#ifndef VERILATED_VATTENTIONSYSTOLICARRAY_H_
#define VERILATED_VATTENTIONSYSTOLICARRAY_H_  // guard

#include "verilated.h"

class VattentionSystolicArray__Syms;
class VattentionSystolicArray___024root;
class VerilatedVcdC;
class VattentionSystolicArray_systolic_array__N4_M4;


// This class is the main interface to the Verilated model
class alignas(VL_CACHE_LINE_BYTES) VattentionSystolicArray VL_NOT_FINAL : public VerilatedModel {
  private:
    // Symbol table holding complete model state (owned by this class)
    VattentionSystolicArray__Syms* const vlSymsp;

  public:

    // CONSTEXPR CAPABILITIES
    // Verilated with --trace?
    static constexpr bool traceCapable = true;

    // PORTS
    // The application code writes and reads these signals to
    // propagate new values into/out from the Verilated model.
    VL_IN8(&clk,0,0);
    VL_IN8(&reset,0,0);
    VL_IN8(&valid_input,0,0);
    VL_OUT8(&valid_result,0,0);
    double &softmax_temp;
    VlUnpacked<VlUnpacked<double, 4>, 4> &K_matrix;
    VlUnpacked<double, 33> &factorial_arr;
    VlUnpacked<VlUnpacked<double, 4>, 4> &Q_matrix;
    VlUnpacked<VlUnpacked<double, 4>, 4> &V_matrix;
    VlUnpacked<VlUnpacked<double, 4>, 4> &Q_mult_K;
    VlUnpacked<VlUnpacked<double, 4>, 4> &exp_Q_mult_K;
    VlUnpacked<VlUnpacked<double, 4>, 4> &exp_K_mult_Q_mult_V;
    VlUnpacked<VlUnpacked<double, 4>, 4> &attention;

    // CELLS
    // Public to allow access to /* verilator public */ items.
    // Otherwise the application code can consider these internals.
    VattentionSystolicArray_systolic_array__N4_M4* const __PVT__attentionSystolicArray__DOT__K_mult_Q_systolicArray;
    VattentionSystolicArray_systolic_array__N4_M4* const __PVT__attentionSystolicArray__DOT__V_mult_systolicArray;

    // Root instance pointer to allow access to model internals,
    // including inlined /* verilator public_flat_* */ items.
    VattentionSystolicArray___024root* const rootp;

    // CONSTRUCTORS
    /// Construct the model; called by application code
    /// If contextp is null, then the model will use the default global context
    /// If name is "", then makes a wrapper with a
    /// single model invisible with respect to DPI scope names.
    explicit VattentionSystolicArray(VerilatedContext* contextp, const char* name = "TOP");
    explicit VattentionSystolicArray(const char* name = "TOP");
    /// Destroy the model; called (often implicitly) by application code
    virtual ~VattentionSystolicArray();
  private:
    VL_UNCOPYABLE(VattentionSystolicArray);  ///< Copying not allowed

  public:
    // API METHODS
    /// Evaluate the model.  Application must call when inputs change.
    void eval() { eval_step(); }
    /// Evaluate when calling multiple units/models per time step.
    void eval_step();
    /// Evaluate at end of a timestep for tracing, when using eval_step().
    /// Application must call after all eval() and before time changes.
    void eval_end_step() {}
    /// Simulation complete, run final blocks.  Application must call on completion.
    void final();
    /// Are there scheduled events to handle?
    bool eventsPending();
    /// Returns time at next time slot. Aborts if !eventsPending()
    uint64_t nextTimeSlot();
    /// Trace signals in the model; called by application code
    void trace(VerilatedTraceBaseC* tfp, int levels, int options = 0) { contextp()->trace(tfp, levels, options); }
    /// Retrieve name of this model instance (as passed to constructor).
    const char* name() const;

    // Abstract methods from VerilatedModel
    const char* hierName() const override final;
    const char* modelName() const override final;
    unsigned threads() const override final;
    /// Prepare for cloning the model at the process level (e.g. fork in Linux)
    /// Release necessary resources. Called before cloning.
    void prepareClone() const;
    /// Re-init after cloning the model at the process level (e.g. fork in Linux)
    /// Re-allocate necessary resources. Called after cloning.
    void atClone() const;
    std::unique_ptr<VerilatedTraceConfig> traceConfig() const override final;
  private:
    // Internal functions - trace registration
    void traceBaseModel(VerilatedTraceBaseC* tfp, int levels, int options);
};

#endif  // guard
