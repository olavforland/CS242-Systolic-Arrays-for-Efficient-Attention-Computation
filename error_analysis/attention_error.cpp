// attention_error.cpp

#include "../common/common.h"
#include "params.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <numeric> // For std::accumulate
#include <iomanip> // For formatting the output
#include <sys/stat.h> // For checking file existence

#define MAX_SIM_TIME 1000
#define RESET_NEG_EDGE 5
#define VERIF_START_TIME 7

// Matrices
std::vector<std::vector<float>> matrix_Q;
std::vector<std::vector<float>> matrix_K;
std::vector<std::vector<float>> matrix_V;
std::vector<std::vector<float>> matrix_attention;

vluint64_t sim_time = 0;

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


// Function to log aggregated MAPE statistics
void verifyAndLogError(float mean_mape, float std_dev_mape, int k_val,
                       float magnitude_scale_val) {

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
        outfile << "K,MagnitudeScale,MeanMAPE,StdDevMAPE\n";
    }

    // Write aggregated results
    outfile << k_val << "," << magnitude_scale_val << "," << mean_mape << "," << std_dev_mape << "\n";
    outfile.close();

    // Print results to console
    std::cout << "\nK: " << k_val << ", Magnitude Scale: " << magnitude_scale_val
              << ", Mean MAPE: " << mean_mape << "%, Standard Deviation of MAPE: " << std_dev_mape << "%\n";
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
    // case 'Q':
    //     displayMatrix("Matrix Q", N, d, [](int i, int j) { return matrix_Q[i][j]; });
    //     break;
    // case 'K':
    //     displayMatrix("Matrix K", N, d, [](int i, int j) { return matrix_K[i][j]; });
    //     break;
    // case 'V':
    //     displayMatrix("Matrix V", N, d, [](int i, int j) { return matrix_V[i][j]; });
    //     break;
    // case 'C':
    //     displayMatrix("Expected QK^T Matrix", N, N, [](int i, int j) { return matrix_Q_mult_K[i][j]; });
    //     break;
    // case 'R':
    //     displayMatrix("Received QK^T Matrix", N, N, [dut](int i, int j) { return dut->Q_mult_K[i][j]; });
    //     break;
    // case 'D':
    //     displayMatrix("Expected exp(QK^T) Matrix", N, N, [](int i, int j) { return matrix_Q_mult_K_exp[i][j]; });
    //     break;
    // case 'E':
    //     displayMatrix("Received exp(QK^T) Matrix", N, N, [dut](int i, int j) { return dut->exp_Q_mult_K[i][j]; });
    //     break;
    // case 'G':
    //     displayMatrix("Expected exp(QK^T)V Matrix", N, d, [](int i, int j) { return matrix_Q_mult_K_exp_mult_V[i][j]; });
    //     break;
    // case 'H':
    //     displayMatrix("Received exp(QK^T)V Matrix", N, d, [dut](int i, int j) { return dut->exp_K_mult_Q_mult_V[i][j]; });
    //     break;
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
        bool incorrect = false;
        calculateResultMatrix(N, d, matrix_Q, matrix_K, matrix_V, matrix_attention);
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

int main(int argc, char** argv, char** env) {
    
    srand(42); // Fixed seed for reproducibility

    // Number of runs for error analysis
    int num_runs = 100;

    

    // Vector to store MAPE values
    std::vector<float> mape_results;

    for (int run = 0; run < num_runs; ++run) {
        Verilated::commandArgs(argc, argv);
        VattentionSystolicArray* dut = new VattentionSystolicArray;


        // Simulation parameters
        const float maxValue = 1.0f;
        const float magnitude_scale_val = MAGNITUDE_SCALE; // From params.h
        const int assertValidInput = (3 * N) + 3;
        const int reset_neg_edge = RESET_NEG_EDGE;
        sim_time = 0;
        vluint64_t posedge_cnt = 0;

        

        // Initialize matrices
        initializeInputMatrices(N, d, maxValue, magnitude_scale_val, matrix_Q, matrix_K, matrix_V);
        loadWeights(dut, N, d, matrix_K, matrix_V);
        calculateResultMatrix(N, d, matrix_Q, matrix_K, matrix_V, matrix_attention);

        // Initialize factorial array and softmax_temp
        create_factorial_arr(K, dut);
        dut->softmax_temp = (float)1.0 / std::sqrt(d);

        // Simulation loop
        while (sim_time < MAX_SIM_TIME) {
            dut_reset(dut, sim_time);
            dut->clk ^= 1; // Toggle clock

            dut->eval(); // Evaluate DUT

            if (dut->clk == 1) {
                ++posedge_cnt;
                toggle_validInput(dut, sim_time, posedge_cnt, assertValidInput, reset_neg_edge);
                driveInputMatrices(dut, sim_time, posedge_cnt, assertValidInput, reset_neg_edge, N, d, matrix_Q);
                // verifyOutputMatrix(dut);

                // verifyAndLogError(dut, N, d, matrix_attention, magnitude_scale_val, K, sim_time);

                // Collect MAPE 
                if (dut->valid_result == 1 && sim_time >= VERIF_START_TIME) {
                    std::vector<std::vector<float>> calculated_attention;
                    float mape = computeMAPE(N, d, matrix_attention, dut, calculated_attention);
                    mape_results.push_back(mape);

            }

            ++sim_time;
            }
        }

        // Clean up
        delete dut;


    }


    std::sort(mape_results.begin(), mape_results.end());

    // Remove the 10 smallest and 10 largest values
    mape_results.erase(mape_results.begin(), mape_results.begin() + 10); // Remove smallest 10
    mape_results.erase(mape_results.end() - 10, mape_results.end());     // Remove largest 10

    // Update the number of runs after trimming
    num_runs = mape_results.size();

    // Calculate mean and standard deviation of MAPE
    float mean_mape = std::accumulate(mape_results.begin(), mape_results.end(), 0.0f) / num_runs;

    float variance = 0.0f;
    for (float mape : mape_results) {
        variance += (mape - mean_mape) * (mape - mean_mape);
    }
    float std_dev_mape = std::sqrt(variance / num_runs);

    // Log
    verifyAndLogError(mean_mape, std_dev_mape, K, MAGNITUDE_SCALE);


    return 0;
}
