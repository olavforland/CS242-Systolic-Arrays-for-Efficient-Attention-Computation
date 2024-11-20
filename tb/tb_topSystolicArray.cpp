#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <vector>

#include "VtopSystolicArray.h" // Verilated DUT.
#include <verilated.h>         // Common verilator routines.
#include <verilated_vcd_c.h>   // Write waverforms to a VCD file.



#define MAX_SIM_TIME 1000 // Number of clk edges.
#define RESET_NEG_EDGE 5  // Clk edge number to deassert arst.
#define VERIF_START_TIME 7

#include "params.h" // Include the automatically generated header


// Max value of an element.
const float maxValue = 1.0f; //static_cast<float>(std::pow(2, 8));
// Number of clk cycles before validInput should be asserted.
const int assertValidInput = (3 * N) + 3;

vluint64_t sim_time = 0;
vluint64_t posedge_cnt = 0;

// Factorial array
// int factorial_arr[K + 1];


float matrix_Q[N][N];
float matrix_K[N][N];
float matrix_V[N][N];
float matrix_Q_mult_K[N][N];
float matrix_Q_mult_K_exp[N][N];
float matrix_Q_mult_K_exp_mult_V[N][N];

// Assert arst only on the first clock edge.
// Note: By default all signals are initialized to 0, so there's no need to
// drive the other inputs to '0.
void dut_reset(VtopSystolicArray *dut) {
  dut->reset = 0;

  if ((sim_time > 2) && (sim_time < RESET_NEG_EDGE)) {
    dut->reset = 1;
  }
}

// Assert validInput after every 10 clk cycles (after reset).
void toggle_validInput(VtopSystolicArray *dut) {
  dut->valid_input = 0;

  if ((posedge_cnt % assertValidInput == 0) && (sim_time >= RESET_NEG_EDGE)) {
    dut->valid_input = 1;
  }
}

void displayMatrix(char matrix, VtopSystolicArray *dut) {
  
  if (matrix == 'Q') {
    std::cout << std::endl;
    std::cout << "Matrix Q " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << std::setprecision(3) << matrix_Q[i][j] << "\t";
      }
      std::cout << std::endl;
    }
  } else if (matrix == 'K') {
    std::cout << std::endl;
    std::cout << "Matrix K " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << std::setprecision(3) << matrix_K[i][j] << "\t";
      }
      std::cout << std::endl;
    }
  } else if (matrix == 'V') {
    std::cout << std::endl;
    std::cout << "Matrix V " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << std::setprecision(3) << matrix_V[i][j] << "\t";
      }
      std::cout << std::endl;
    }
  } else if (matrix == 'C') {
    std::cout << std::endl;
    std::cout << "Expected QK^T Matrix " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << std::setprecision(3) <<  matrix_Q_mult_K[i][j] << "\t";
      }
      std::cout << std::endl;
    }
  } else if (matrix == 'R') {
    std::cout << std::endl;
    std::cout << "Received QK^T Matrix " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << std::setprecision(3) << dut->K_mult_Q_out[i][j]
                  << "\t";
      }
      std::cout << std::endl;
    }
  }
  else if (matrix == 'E') {
    std::cout << std::endl;
    std::cout << "Received exp(QK^T) Matrix " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << dut->exponentiation_out[i][j]
                  << "\t";
      }
      std::cout << std::endl;
    }
  }
  else if (matrix == 'D') {
    std::cout << std::endl;
    std::cout << "Expected exp(QK^T) Matrix " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << matrix_Q_mult_K_exp[i][j]
                  << "\t";
      }
      std::cout << std::endl;
    }
  }
  else if (matrix == 'G') {
    std::cout << std::endl;
    std::cout << "Expected exp(QK^T)V Matrix " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << matrix_Q_mult_K_exp_mult_V[i][j]
                  << "\t";
      }
      std::cout << std::endl;
    }
  }
  else if (matrix == 'H') {
    std::cout << std::endl;
    std::cout << "Received exp(QK^T)V Matrix " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << dut->exp_K_mult_Q_mult_V[i][j]
                  << "\t";
      }
      std::cout << std::endl;
    }
  }
  
}

void create_factorial_arr(int K_val, VtopSystolicArray *dut) {
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
    for (int j = 0; j < N; j++) {
      matrix_Q[i][j] = static_cast<float>(rand()) / RAND_MAX * maxValue;

      // matrix_Q[i][j] = static_cast<float>(rand() % maxValue);
      // std::cout << matrix_Q[i][j] << " ";
      matrix_K[i][j] = static_cast<float>(rand()) / RAND_MAX * maxValue;
      // matrix_K[i][j] = rand() % maxValue * 1.0;

      matrix_V[i][j] = static_cast<float>(rand()) / RAND_MAX * maxValue;
    }
  }
  displayMatrix('Q', nullptr);
  displayMatrix('K', nullptr);
  displayMatrix('V', nullptr);
}
void loadWeights(VtopSystolicArray *dut) {
    // Provide weight values to 'weight_in'
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            dut->K_matrix[i][j] = matrix_K[i][j];
            dut->V_matrix[i][j] = matrix_V[i][j];
        }
    }
}


void driveInputMatrices(VtopSystolicArray *dut) {
    // Persistent input matrix storage
    static float data[N][N] = {0}; // Keeps track of the current state of inputs

    if (sim_time < RESET_NEG_EDGE) {
        // During reset, ensure data_in is zeroed
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                dut->Q_matrix[i][j] = 0;
            }
        }
        return; // Skip the rest of the logic during reset
    }

    // After reset, persist input matrix values
    if (posedge_cnt % assertValidInput == 0) {
        // Assign `matrix_Q` to the persistent data buffer
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                data[i][j] = matrix_Q[i][j];
                dut->Q_matrix[i][j] = data[i][j];
            }
        }
    }
}


void calculateResultMatrix() {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      matrix_Q_mult_K[i][j] = 0.0;

      for (int k = 0; k < N; ++k) {
        matrix_Q_mult_K[i][j] += matrix_Q[i][k] * matrix_K[j][k];
      }
      matrix_Q_mult_K_exp[i][j] = std::exp(matrix_Q_mult_K[i][j]);
    }
  }
    // Step 2: Compute the matrix product of exp(QK^T) with V
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      matrix_Q_mult_K_exp_mult_V[i][j] = 0.0; // Initialize result matrix element

      for (int k = 0; k < N; ++k) {
        matrix_Q_mult_K_exp_mult_V[i][j] += matrix_Q_mult_K_exp[i][k] * matrix_V[k][j];
      }
    }
  }
}

void verifyOutputMatrix(VtopSystolicArray *dut) {
  if ((dut->valid_result == 1) && (sim_time >= VERIF_START_TIME)) {
    calculateResultMatrix();
    displayMatrix('C', dut);

    // Note: Verilator represents the output matrix as n^2 bit arrays.
    bool incorrect = false;

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (std::abs(dut->K_mult_Q_out[i][j] - matrix_Q_mult_K[i][j]) > 1e-2) {
            incorrect = true;
          }
      }
    }
    displayMatrix('R', dut);
    displayMatrix('E', dut);
    displayMatrix('D', dut);

    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (std::abs(dut->exp_K_mult_Q_mult_V[i][j] - matrix_Q_mult_K_exp_mult_V[i][j]) > 1e-2) {
            incorrect = true;
          }
      }
    }
    displayMatrix('G', dut);
    displayMatrix('H', dut);
    if (incorrect) {
      std::cout << std::endl;
      std::cerr << "ERROR: output matrix received is incorrect." << std::endl;
      displayMatrix('R', dut);
      std::cout << " simtime: " << (int)sim_time << std::endl;
      std::cout << "*******************************************" << std::endl;
      // exit(EXIT_FAILURE);
    }
  }
}


int main(int argc, char **argv, char **env) {
  srand(time(NULL));
  Verilated::commandArgs(argc, argv);
  VtopSystolicArray *dut = new VtopSystolicArray; // Instantiate DUT.

  // Initialize matrices
  initializeInputMatrices();
  create_factorial_arr(K, dut);
  calculateResultMatrix(); // Precompute the expected result
  loadWeights(dut);

  // {{{ Set-up waveform dumping.

  Verilated::traceEverOn(true);
  VerilatedVcdC *m_trace = new VerilatedVcdC;
  dut->trace(m_trace, 5);
  m_trace->open("waveform.vcd");

  // }}} Set-up waveform dumping

  while (sim_time < MAX_SIM_TIME) {
    dut_reset(dut);

    dut->clk ^= 1; // Toggle clk to create pos and neg edge.

    dut->eval(); // Evaluate all the signals in the DUT on each clock edge.

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
}
