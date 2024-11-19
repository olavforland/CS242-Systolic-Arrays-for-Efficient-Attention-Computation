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
const int maxValue = std::pow(2, 8);
// Number of clk cycles before validInput should be asserted.
const int assertValidInput = (3 * N) + 3;

vluint64_t sim_time = 0;
vluint64_t posedge_cnt = 0;

uint8_t matrixA[N][N];
uint8_t matrixB[N][N];
uint32_t matrixC[N][N];

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
  if (matrix == 'A') {
    std::cout << std::endl;
    std::cout << "Matrix A " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << std::hex << static_cast<int>(matrixA[i][j]) << "\t";
      }
      std::cout << std::endl;
    }
  } else if (matrix == 'B') {
    std::cout << std::endl;
    std::cout << "Matrix B " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << std::hex << static_cast<int>(matrixB[i][j]) << "\t";
      }
      std::cout << std::endl;
    }
  } else if (matrix == 'C') {
    std::cout << std::endl;
    std::cout << "Expected result Matrix:" << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << std::hex << static_cast<int>(matrixC[i][j]) << "\t";
      }
      std::cout << std::endl;
    }
  } else if (matrix == 'R') {
    std::cout << std::endl;
    std::cout << "Received matrix " << std::endl;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        std::cout << std::hex << static_cast<int>(dut->result_out[i][j])
                  << "\t";
      }
      std::cout << std::endl;
    }
  }
}


void initializeInputMatrices() {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      matrixA[i][j] = rand() % maxValue;
      matrixB[i][j] = rand() % maxValue;
    }
  }
  displayMatrix('A', nullptr);
  displayMatrix('B', nullptr);
}
void loadWeights(VtopSystolicArray *dut) {
    // Provide weight values to 'weight_in'
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            dut->weight_in[j][i] = matrixB[i][j];
        }
    }
}


void driveInputMatrices(VtopSystolicArray *dut) {
    // Persistent input matrix storage
    static int data[N][N] = {0}; // Keeps track of the current state of inputs

    if (sim_time < RESET_NEG_EDGE) {
        // During reset, ensure data_in is zeroed
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                dut->data_in[i][j] = 0;
            }
        }
        return; // Skip the rest of the logic during reset
    }

    // After reset, persist input matrix values
    if (posedge_cnt % assertValidInput == 0) {
        // Assign `matrixA` to the persistent data buffer
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                data[i][j] = matrixA[i][j];
                dut->data_in[i][j] = data[i][j];
            }
        }
    }
}


void calculateResultMatrix() {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      matrixC[i][j] = 0;

      for (int k = 0; k < N; ++k) {
        matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
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
        //if (dut->result_out[(N * i) + j] != matrixC[i][j]) {
        if (dut->result_out[i][j] != matrixC[i][j]) {
          incorrect = true;
        }
      }
    }
    displayMatrix('R', dut);
    if (incorrect) {
      std::cout << std::endl;
      std::cerr << "ERROR: output matrix received is incorrect." << std::endl;
      displayMatrix('R', dut);
      std::cout << " simtime: " << (int)sim_time << std::endl;
      std::cout << "*******************************************" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}


int main(int argc, char **argv, char **env) {
  srand(time(NULL));
  Verilated::commandArgs(argc, argv);
  VtopSystolicArray *dut = new VtopSystolicArray; // Instantiate DUT.

  // Initialize matrices
  initializeInputMatrices();
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
