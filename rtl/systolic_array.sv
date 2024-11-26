`timescale 1ns / 1ps
`default_nettype none

module systolic_array #(
    parameter int unsigned N, // Outer dimension of weights
    parameter int unsigned M  // Inner dimension that joins the matrices
)(
    input logic clk,
    input logic reset,
    input logic weight_load_enable,
    input logic doProcess,
    // Weight input for each PE
    input real weight_in [0:N-1][0:M-1],
    // Input for the first row of PEs
    input real data_in [0:M-1],
    // Resulting output matrix
    output real result_out [0:N-1]
);
    // Wires for data flowing upwards
    real data_wire [0:N][0:M-1];
    
    // Wires for partial sums flowing to the right
    real partial_sum_wire [0:N-1][0:M];

    // Boundary conditions for data flow from the input to the first row
    generate
        for (genvar j = 0; j < M; j++) begin : data_flow
            assign data_wire[N][j] = data_in[j];
        end
    endgenerate

    // Boundary conditions for partial sums starting at zero for the first column
    generate
        for (genvar i = 0; i < N; i++) begin : sum_flow
            assign partial_sum_wire[i][0] = '0;
        end
    endgenerate

    // Instantiate the PEs and connect inter-PE signals
    generate
        for (genvar i = 0; i < N; i++) begin : row
            for (genvar j = 0; j < M; j++) begin : col
                pe #(
                ) pe (
                    .clk(clk),
                    .reset(reset),
                    .doProcess(doProcess),
                    .weight_load_enable(weight_load_enable),
                    .weight_in(weight_in[i][j]),
                    .data_in(data_wire[i+1][j]),
                    .partial_sum_in(partial_sum_wire[i][j]),
                    .data_out(data_wire[i][j]),
                    .partial_sum_out(partial_sum_wire[i][j+1])
                );
            end
        end
    endgenerate

    generate
        for (genvar i = 0; i < N; i++) begin : result
            assign result_out[i] = partial_sum_wire[i][M];
        end
    endgenerate

endmodule

`default_nettype wire
