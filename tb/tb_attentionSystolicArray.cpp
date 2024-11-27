#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

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

float matrix_Q[N][d];
float matrix_K[N][d];
float matrix_V[N][d];
float matrix_Q_mult_K[N][N];
float matrix_Q_mult_K_exp[N][N];
float matrix_Q_mult_K_exp_mult_V[N][d];
float matrix_attention_norm[N];
float matrix_attention[N][d];

// Assert arst only on the first clock edge.
void dut_reset(VattentionSystolicArray *dut) {
    dut->reset = 0;
    if ((sim_time > 2) && (sim_time < RESET_NEG_EDGE)) {
        dut->reset = 1;
    }
}

// Assert validInput after every 'assertValidInput' clk cycles (after reset).
void toggle_validInput(VattentionSystolicArray *dut) {
    dut->valid_input = 0;
    if ((posedge_cnt % assertValidInput == 0) && (sim_time >= RESET_NEG_EDGE)) {
        dut->valid_input = 1;
    }
}

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

void create_factorial_arr(int K_val, VattentionSystolicArray *dut) {
    int fact = 1;
    dut->factorial_arr[0] = 1;
    for (int i = 1; i <= K; i++) {
        fact = fact * i;
        dut->factorial_arr[i] = fact;
    }
}

void initializeInputMatrices() {
    srand(42);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < d; j++) {
            matrix_Q[i][j] = static_cast<float>(rand()) / RAND_MAX * maxValue;
            matrix_K[i][j] = static_cast<float>(rand()) / RAND_MAX * maxValue;
            matrix_V[i][j] = static_cast<float>(rand()) / RAND_MAX * maxValue;
        }
    }
    displayMatrix('Q', nullptr);
    displayMatrix('K', nullptr);
    displayMatrix('V', nullptr);
}

void loadWeights(VattentionSystolicArray *dut) {
    // Provide weight values to 'weight_in'
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < d; ++j) {
            dut->K_matrix[i][j] = matrix_K[i][j];
            dut->V_matrix[i][j] = matrix_V[i][j];
        }
    }
}

void driveInputMatrices(VattentionSystolicArray *dut) {
    if (sim_time < RESET_NEG_EDGE) {
        // During reset, ensure Q_matrix is zeroed
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < d; j++) {
                dut->Q_matrix[i][j] = 0;
            }
        }
    } else if (posedge_cnt % assertValidInput == 0) {
        // Assign matrix_Q to dut->Q_matrix
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < d; j++) {
                dut->Q_matrix[i][j] = matrix_Q[i][j];
            }
        }
    }
}

void calculateResultMatrix() {
    for (int i = 0; i < N; i++) {
        matrix_attention_norm[i] = 0.0;
        for (int j = 0; j < N; j++) {
            matrix_Q_mult_K[i][j] = 0.0;
            for (int k = 0; k < d; ++k) {
                matrix_Q_mult_K[i][j] += (matrix_Q[i][k] * matrix_K[j][k]);
            }
            matrix_Q_mult_K_exp[i][j] = std::exp(matrix_Q_mult_K[i][j] / std::sqrt((float) d));
            matrix_attention_norm[i] += matrix_Q_mult_K_exp[i][j];
        }
    }
    // Compute the matrix product of exp(QK^T) with V
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < d; j++) {
            matrix_Q_mult_K_exp_mult_V[i][j] = 0.0;
            for (int k = 0; k < N; ++k) {
                matrix_Q_mult_K_exp_mult_V[i][j] += (matrix_Q_mult_K_exp[i][k] * matrix_V[k][j]);
            }
            matrix_attention[i][j] = matrix_Q_mult_K_exp_mult_V[i][j] / matrix_attention_norm[i];
        }
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
        calculateResultMatrix();
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
    initializeInputMatrices();
    create_factorial_arr(K, dut);
    calculateResultMatrix(); // Precompute the expected result
    loadWeights(dut);
    dut->softmax_temp = (float) 1.0 / std::sqrt(d);

    // Set-up waveform dumping.
    Verilated::traceEverOn(true);
    VerilatedVcdC *m_trace = new VerilatedVcdC;
    dut->trace(m_trace, 5);
    m_trace->open("waveform.vcd");

    while (sim_time < MAX_SIM_TIME) {
        dut_reset(dut);
        dut->clk ^= 1; // Toggle clk to create pos and neg edge.
        dut->eval();   // Evaluate all the signals in the DUT on each clock edge.

        if (dut->clk == 1) {
            posedge_cnt++;
            toggle_validInput(dut);
            driveInputMatrices(dut);
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
