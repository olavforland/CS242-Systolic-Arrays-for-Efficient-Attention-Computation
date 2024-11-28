#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../common/common.h"

#include "VattentionSystolicArray.h" // Verilated DUT.
#include <verilated.h>               // Common verilator routines.
#include <verilated_vcd_c.h>         // Write waverforms to a VCD file.

#define MAX_SIM_TIME 1000 // Number of clk edges.
#define RESET_NEG_EDGE 5  // Clk edge number to deassert arst.
#define VERIF_START_TIME 7

#include "params.h" // Include the automatically generated header

// Max value of an element.
const float maxValue = 1.0f;
// Number of clk cycles before validInput should be asserted.
const int assertValidInput = (3 * N) + 3;

vluint64_t sim_time = 0;
vluint64_t posedge_cnt = 0;

// float matrix_Q[N][d];
// float matrix_K[N][d];
// float matrix_V[N][d];
float matrix_Q_mult_K[N][N];
float matrix_Q_mult_K_exp[N][N];
float matrix_Q_mult_K_exp_mult_V[N][d];
float matrix_attention_norm[N];
// float matrix_attention[N][d];

std::vector<std::vector<float>> matrix_Q;
std::vector<std::vector<float>> matrix_K;
std::vector<std::vector<float>> matrix_V;
std::vector<std::vector<float>> matrix_attention;


void displayMatrix(const std::string &name, int rows, int cols, std::function<float(int, int)> getData) {
    std::cout << "\n" << name << std::endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << std::setprecision(3) << getData(i, j) << "\t";
        }
        std::cout << std::endl;
    }
}

void displayMatrix(char matrix, VattentionSystolicArray *dut) {
    switch (matrix) {
    case 'Q':
        displayMatrix("Matrix Q", N, d, [](int i, int j) { return matrix_Q[i][j]; });
        break;
    case 'K':
        displayMatrix("Matrix K", N, d, [](int i, int j) { return matrix_K[i][j]; });
        break;
    case 'V':
        displayMatrix("Matrix V", N, d, [](int i, int j) { return matrix_V[i][j]; });
        break;
    case 'C':
        displayMatrix("Expected QK^T Matrix", N, N, [](int i, int j) { return matrix_Q_mult_K[i][j]; });
        break;
    case 'R':
        displayMatrix("Received QK^T Matrix", N, N, [dut](int i, int j) { return dut->Q_mult_K[i][j]; });
        break;
    case 'D':
        displayMatrix("Expected exp(QK^T) Matrix", N, N, [](int i, int j) { return matrix_Q_mult_K_exp[i][j]; });
        break;
    case 'E':
        displayMatrix("Received exp(QK^T) Matrix", N, N, [dut](int i, int j) { return dut->exp_Q_mult_K[i][j]; });
        break;
    case 'G':
        displayMatrix("Expected exp(QK^T)V Matrix", N, d, [](int i, int j) { return matrix_Q_mult_K_exp_mult_V[i][j]; });
        break;
    case 'H':
        displayMatrix("Received exp(QK^T)V Matrix", N, d, [dut](int i, int j) { return dut->exp_K_mult_Q_mult_V[i][j]; });
        break;
    case 'I':
        displayMatrix("Expected Attention Matrix", N, d, [](int i, int j) { return matrix_attention[i][j]; });
        break;
    case 'J':
        displayMatrix("Received Attention Matrix", N, d, [dut](int i, int j) { return dut->attention[i][j]; });
        break;
    default:
        std::cerr << "Unknown matrix code: " << matrix << std::endl;
    }
}

void compareMatrices(int rows, int cols, std::function<float(int, int)> data1,
                     std::function<float(int, int)> data2, bool &incorrect) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (std::abs(data1(i, j) - data2(i, j)) > 1e-2) {
                incorrect = true;
            }
        }
    }
}

void verifyOutputMatrix(VattentionSystolicArray *dut) {
    if ((dut->valid_result == 1) && (sim_time >= VERIF_START_TIME)) {
        calculateResultMatrix(N, d, matrix_Q, matrix_K, matrix_V, matrix_attention);
        displayMatrix('C', dut);

        bool incorrect = false;

        compareMatrices(N, N, [](int i, int j) { return matrix_Q_mult_K[i][j]; },
                        [dut](int i, int j) { return dut->Q_mult_K[i][j]; }, incorrect);
        displayMatrix('R', dut);
        displayMatrix('D', dut);
        displayMatrix('E', dut);

        compareMatrices(N, d, [](int i, int j) { return matrix_Q_mult_K_exp_mult_V[i][j]; },
                        [dut](int i, int j) { return dut->exp_K_mult_Q_mult_V[i][j]; }, incorrect);
        displayMatrix('G', dut);
        displayMatrix('H', dut);

        compareMatrices(N, d, [](int i, int j) { return matrix_attention[i][j]; },
                        [dut](int i, int j) { return dut->attention[i][j]; }, incorrect);
        displayMatrix('I', dut);
        displayMatrix('J', dut);

        if (incorrect) {
            std::cerr << "\nERROR: output matrix received is incorrect." << std::endl;
            std::cout << " simtime: " << (int)sim_time << std::endl;
            std::cout << "*******************************************" << std::endl;
        }
    }
}

int main(int argc, char **argv, char **env) {
    Verilated::commandArgs(argc, argv);
    VattentionSystolicArray *dut = new VattentionSystolicArray; // Instantiate DUT.

    // Initialize matrices
    initializeInputMatrices(N, d, maxValue, 1.0, matrix_Q, matrix_K, matrix_V);
    create_factorial_arr(K, dut);
    calculateResultMatrix(N, d, matrix_Q, matrix_K, matrix_V, matrix_attention);

    loadWeights(dut, N, d, matrix_K, matrix_V);
    dut->softmax_temp = (float) 1.0 / std::sqrt(d);

    // Set-up waveform dumping.
    Verilated::traceEverOn(true);
    VerilatedVcdC *m_trace = new VerilatedVcdC;
    dut->trace(m_trace, 5);
    m_trace->open("waveform.vcd");

    while (sim_time < MAX_SIM_TIME) {
        dut_reset(dut, sim_time);
        dut->clk ^= 1; // Toggle clk to create pos and neg edge.
        dut->eval();   // Evaluate all the signals in the DUT on each clock edge.

        if (dut->clk == 1) {
            posedge_cnt++;
            toggle_validInput(dut, sim_time, posedge_cnt, assertValidInput, RESET_NEG_EDGE);
            driveInputMatrices(dut, sim_time, posedge_cnt, assertValidInput, RESET_NEG_EDGE, N, d, matrix_Q);
            verifyOutputMatrix(dut);
        }

        // Write all the traced signal values into the waveform dump file.
        m_trace->dump(sim_time);
        sim_time++;
    }
    m_trace->close();
    delete dut;
    delete m_trace;
    return 0;
}
