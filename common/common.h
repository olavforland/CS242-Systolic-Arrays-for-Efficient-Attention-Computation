// common.h
#ifndef COMMON_H
#define COMMON_H

#include "VattentionSystolicArray.h"
#include <verilated.h>
#include <verilated_vcd_c.h>

// Function declarations

void dut_reset(VattentionSystolicArray* dut, vluint64_t sim_time);
void toggle_validInput(VattentionSystolicArray* dut, vluint64_t sim_time, vluint64_t posedge_cnt, int assertValidInput, int reset_neg_edge);
void initializeInputMatrices(int N, int d, float maxValue, float magnitude_scale,
                             std::vector<std::vector<float>>& matrix_Q,
                             std::vector<std::vector<float>>& matrix_K,
                             std::vector<std::vector<float>>& matrix_V);
void loadWeights(VattentionSystolicArray* dut, int N, int d,
                 const std::vector<std::vector<float>>& matrix_K,
                 const std::vector<std::vector<float>>& matrix_V);
void driveInputMatrices(VattentionSystolicArray* dut, vluint64_t sim_time, vluint64_t posedge_cnt,
                        int assertValidInput, int reset_neg_edge,
                        int N, int d, const std::vector<std::vector<float>>& matrix_Q);
void calculateResultMatrix(int N, int d,
                           const std::vector<std::vector<float>>& matrix_Q,
                           const std::vector<std::vector<float>>& matrix_K,
                           const std::vector<std::vector<float>>& matrix_V,
                           std::vector<std::vector<float>>& matrix_attention,
                           std::vector<std::vector<float>>& matrix_Q_mult_K = *(new std::vector<std::vector<float>>()),
                           std::vector<std::vector<float>>& matrix_Q_mult_K_exp = *(new std::vector<std::vector<float>>()),
                           std::vector<std::vector<float>>& matrix_Q_mult_K_exp_mult_V = *(new std::vector<std::vector<float>>()),
                           std::vector<float>& matrix_attention_norm = *(new std::vector<float>()));
void create_factorial_arr(int k_val, VattentionSystolicArray* dut);

#endif // COMMON_H