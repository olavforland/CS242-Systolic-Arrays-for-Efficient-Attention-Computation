// attention_error.cpp

#include "../common/common.h"
#include "params.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip> // For formatting the output
#include <sys/stat.h> // For checking file existence

#define MAX_SIM_TIME 1000
#define RESET_NEG_EDGE 5
#define VERIF_START_TIME 7

// Function to compute Mean Absolute Percentage Error (MAPE)
float computeMAPE(int n_val, int d_val,
                  const std::vector<std::vector<float>>& expected_attention,
                  VattentionSystolicArray* dut,
                  std::vector<std::vector<float>>& calculated_attention) {
    float total_percentage_error = 0.0f;
    int count = 0;

    calculated_attention.resize(n_val, std::vector<float>(d_val));

    for (int i = 0; i < n_val; ++i) {
        for (int j = 0; j < d_val; ++j) {
            float expected = expected_attention[i][j];
            float actual = dut->attention[i][j];
            calculated_attention[i][j] = actual;

            if (expected != 0.0f) {
                float percentage_error = std::abs((expected - actual) / expected) * 100.0f;
                total_percentage_error += percentage_error;
                ++count;
            }
        }
    }

    return (count > 0) ? (total_percentage_error / count) : 0.0f;
}

// Function to verify output, print matrices, and write MAPE to CSV
void verifyAndLogError(VattentionSystolicArray* dut, int n_val, int d_val,
                       const std::vector<std::vector<float>>& expected_attention,
                       float magnitude_scale_val, int k_val, vluint64_t& sim_time) {
    if (dut->valid_result == 1 && sim_time >= VERIF_START_TIME) {
        std::vector<std::vector<float>> calculated_attention;
        float mape = computeMAPE(n_val, d_val, expected_attention, dut, calculated_attention);

        // Ensure results directory exists
        system("mkdir -p results");

        // Check if the CSV file already exists
        std::string csv_filename = "results/error_results.csv";
        std::ifstream infile(csv_filename);
        bool file_exists = infile.good();
        infile.close();

        // Open the CSV file in append mode
        std::ofstream outfile(csv_filename, std::ios_base::app);

        // Write headers only if the file did not exist before
        if (!file_exists) {
            outfile << "K,MagnitudeScale,MAPE\n";
        }

        outfile << k_val << "," << magnitude_scale_val << "," << mape << "\n";
        outfile.close();

        // Set precision for floating-point output
        std::cout << std::fixed << std::setprecision(6);

        // Print MAPE to console
        std::cout << "\nK: " << k_val << ", Magnitude Scale: " << magnitude_scale_val
                  << ", MAPE: " << mape << "%" << std::endl;

        // Print Expected Attention Matrix
        std::cout << "\nExpected Attention Matrix:\n";
        for (int i = 0; i < n_val; ++i) {
            for (int j = 0; j < d_val; ++j) {
                std::cout << std::setw(12) << expected_attention[i][j] << " ";
            }
            std::cout << "\n";
        }

        // Print Calculated Attention Matrix
        std::cout << "\nCalculated Attention Matrix from DUT:\n";
        for (int i = 0; i < n_val; ++i) {
            for (int j = 0; j < d_val; ++j) {
                std::cout << std::setw(12) << calculated_attention[i][j] << " ";
            }
            std::cout << "\n";
        }

        // End simulation
        sim_time = MAX_SIM_TIME; // Terminate simulation
    }
}

int main(int argc, char** argv, char** env) {
    Verilated::commandArgs(argc, argv);
    VattentionSystolicArray* dut = new VattentionSystolicArray;

    // Simulation parameters
    const float maxValue = 1.0f;
    const float magnitude_scale_val = MAGNITUDE_SCALE; // From params.h
    const int assertValidInput = (3 * N) + 3;
    const int reset_neg_edge = RESET_NEG_EDGE;
    vluint64_t sim_time = 0;
    vluint64_t posedge_cnt = 0;

    // Matrices
    std::vector<std::vector<float>> matrix_Q;
    std::vector<std::vector<float>> matrix_K;
    std::vector<std::vector<float>> matrix_V;
    std::vector<std::vector<float>> matrix_attention;

    // Initialize matrices
    initializeInputMatrices(N, d, maxValue, magnitude_scale_val, matrix_Q, matrix_K, matrix_V);
    loadWeights(dut, N, d, matrix_K, matrix_V);
    calculateResultMatrix(N, d, matrix_Q, matrix_K, matrix_V, matrix_attention);

    // Initialize factorial array and softmax_temp
    create_factorial_arr(K, dut);
    dut->softmax_temp = (float)1.0 / std::sqrt(d);

    // Set up waveform dumping (if needed)
    Verilated::traceEverOn(true);
    VerilatedVcdC* m_trace = new VerilatedVcdC;
    dut->trace(m_trace, 5);
    m_trace->open("waveform.vcd");

    // Simulation loop
    while (sim_time < MAX_SIM_TIME) {
        dut_reset(dut, sim_time);
        dut->clk ^= 1; // Toggle clock

        dut->eval(); // Evaluate DUT

        if (dut->clk == 1) {
            ++posedge_cnt;
            toggle_validInput(dut, sim_time, posedge_cnt, assertValidInput, reset_neg_edge);
            driveInputMatrices(dut, sim_time, posedge_cnt, assertValidInput, reset_neg_edge, N, d, matrix_Q);

            verifyAndLogError(dut, N, d, matrix_attention, magnitude_scale_val, K, sim_time);
        }

        m_trace->dump(sim_time); // Dump waveforms
        ++sim_time;
    }

    // Clean up
    m_trace->close();
    delete dut;
    delete m_trace;

    return 0;
}
