// common.cpp
#include "common.h"
#include <vector>
#include <cmath>
#include <cstdlib>

void dut_reset(VattentionSystolicArray* dut, vluint64_t sim_time) {
    dut->reset = (sim_time > 2 && sim_time < 5) ? 1 : 0;
}

void toggle_validInput(VattentionSystolicArray* dut, vluint64_t sim_time, vluint64_t posedge_cnt, int assertValidInput, int reset_neg_edge) {
    dut->valid_input = 0;
    if ((posedge_cnt % assertValidInput == 0) && (sim_time >= reset_neg_edge)) {
        dut->valid_input = 1;
    }
}

void create_factorial_arr(int k_val, VattentionSystolicArray* dut) {
    int fact = 1;
    dut->factorial_arr[0] = 1;
    for (int i = 1; i <= k_val; i++) {
        fact = fact * i;
        dut->factorial_arr[i] = fact;
    }
}

void initializeInputMatrices(int N, int d, float maxValue, float magnitude_scale,
                             std::vector<std::vector<float>>& matrix_Q,
                             std::vector<std::vector<float>>& matrix_K,
                             std::vector<std::vector<float>>& matrix_V) {
    srand(42); // Fixed seed for reproducibility
    matrix_Q.resize(N, std::vector<float>(d));
    matrix_K.resize(N, std::vector<float>(d));
    matrix_V.resize(N, std::vector<float>(d));

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < d; ++j) {
            matrix_Q[i][j] = static_cast<float>(rand()) / RAND_MAX * maxValue * magnitude_scale;
            matrix_K[i][j] = static_cast<float>(rand()) / RAND_MAX * maxValue * magnitude_scale;
            matrix_V[i][j] = static_cast<float>(rand()) / RAND_MAX * maxValue * magnitude_scale;
        }
    }
}

void loadWeights(VattentionSystolicArray* dut, int N, int d,
                 const std::vector<std::vector<float>>& matrix_K,
                 const std::vector<std::vector<float>>& matrix_V) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < d; ++j) {
            dut->K_matrix[i][j] = matrix_K[i][j];
            dut->V_matrix[i][j] = matrix_V[i][j];
        }
    }
}

void driveInputMatrices(VattentionSystolicArray* dut, vluint64_t sim_time, vluint64_t posedge_cnt,
                        int assertValidInput, int reset_neg_edge,
                        int N, int d, const std::vector<std::vector<float>>& matrix_Q) {
    if (sim_time < reset_neg_edge) {
        // During reset, ensure Q_matrix is zeroed
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < d; ++j) {
                dut->Q_matrix[i][j] = 0;
            }
        }
    } else if (posedge_cnt % assertValidInput == 0) {
        // Assign matrix_Q to dut->Q_matrix
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < d; ++j) {
                dut->Q_matrix[i][j] = matrix_Q[i][j];
            }
        }
    }
}

void calculateResultMatrix(int N, int d,
                           const std::vector<std::vector<float>>& matrix_Q,
                           const std::vector<std::vector<float>>& matrix_K,
                           const std::vector<std::vector<float>>& matrix_V,
                           std::vector<std::vector<float>>& matrix_attention) {
    std::vector<std::vector<float>> matrix_Q_mult_K(N, std::vector<float>(N));
    std::vector<std::vector<float>> matrix_Q_mult_K_exp(N, std::vector<float>(N));
    std::vector<float> matrix_attention_norm(N, 0.0f);
    matrix_attention.resize(N, std::vector<float>(d));

    // Compute QK^T and exp(QK^T)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < d; ++k) {
                sum += matrix_Q[i][k] * matrix_K[j][k];
            }
            matrix_Q_mult_K[i][j] = sum;
            matrix_Q_mult_K_exp[i][j] = std::exp(sum / std::sqrt((float) d));
            matrix_attention_norm[i] += matrix_Q_mult_K_exp[i][j];
        }
    }

    // Compute exp(QK^T) * V and normalize
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < d; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < N; ++k) {
                sum += matrix_Q_mult_K_exp[i][k] * matrix_V[k][j];
            }
            matrix_attention[i][j] = sum / matrix_attention_norm[i];
        }
    }
}
