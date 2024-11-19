`timescale 1ns / 1ps
`default_nettype none

module systolic_array_exp #(
    parameter int unsigned K = 4,  // Number of Taylor terms (columns)
    parameter int unsigned N = 4   // Number of rows
)(
    input logic clk,
    input logic reset,
    input logic doProcess,
    input real factorial_arr [0:K], // Factorial terms for each column
    input real data_in [0:N-1], // Input column vector
    output real exp_out [0:N-1], // Output exponent approximation
    output real exp_sum_out     // Output row-wise accumulated sum
);

    // Internal wires for systolic array
    real x_wire [0:N-1][0:K];              // Data flows left to right
    real x_pow_wire [0:N-1][0:K];          // x^k flows left to right
    real taylor_approx_wire [0:N-1][0:K];  // Taylor approx flows left to right
    real accumulated_sum_wire [0:N];       // Accumulated sum flows upwards (final column only)

    // Boundary conditions for x_wire
    generate
        for (genvar i = 0; i < N; i++) begin : boundary
            assign x_wire[i][0] = data_in[i]; // Data input to the leftmost column
            assign x_pow_wire[i][0] = 1.0;          // x^0 = 1.0
            assign taylor_approx_wire[i][0] = 1.0;  // Initial Taylor approximation
        end
    endgenerate

    // Boundary condition for accumulated_sum_wire
    assign accumulated_sum_wire[N] = 0.0; // Bottom boundary of accumulation

    // Instantiate PEs for columns 0 to K-2
    for (genvar i = 0; i < N; i++) begin : row
        for (genvar j = 0; j < K-1; j++) begin : col
            pe_exp pe (
                .clk(clk),
                .reset(reset),
                .doProcess(doProcess),
                .k_factorial(factorial_arr[j+1]),                      // Factorial term
                .x_in(x_wire[i][j]),                  // Data input from the left
                .x_pow_in(x_pow_wire[i][j]),             // x^k input
                .taylor_approx_in(taylor_approx_wire[i][j]), // Taylor sum from the left
                .x_out(x_wire[i][j+1]),               // Data output to the right
                .x_pow_out(x_pow_wire[i][j+1]),          // x^(k+1) output to the right
                .taylor_approx_out(taylor_approx_wire[i][j+1]) // Updated Taylor sum
            );
        end
        // Final column: pe_exp_acc for accumulation and output
        pe_exp_acc pe_acc (
            .clk(clk),
            .reset(reset),
            .doProcess(doProcess),
            .exp_in(taylor_approx_wire[i][K-1]),       // Final Taylor sum
            .exp_sum_in(accumulated_sum_wire[i+1]), // Accumulated sum from below
            .exp_out(taylor_approx_wire[i][K]),     // Taylor sum to the right
            .exp_sum_out(accumulated_sum_wire[i])   // Updated accumulated sum
        );
    end

    // Assign outputs from accumulated_sum_wire (row-wise accumulated sums)
    assign exp_sum_out = accumulated_sum_wire[0];

    generate
        for (genvar i = 0; i < N; i++) begin : assign_exp_out
            assign exp_out[i] = taylor_approx_wire[i][K]; // Final Taylor approximation for each row
        end
    endgenerate


endmodule

`default_nettype wire
