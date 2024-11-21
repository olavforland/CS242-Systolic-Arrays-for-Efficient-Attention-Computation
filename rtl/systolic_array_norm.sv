`timescale 1ns / 1ps
`default_nettype none

module systolic_array_norm #(
    parameter int unsigned N = 4   // Number of elements
)(
    input logic clk,
    input logic reset,
    input logic doProcess,
    input real norm_in,                  // Normalizing factor
    input real data_in [0:N-1],          // Input row vector
    output real data_normed_out [0:N-1]  // Row-wise normalized data output
);

    // Internal wires for systolic array
    real norm_wire [0:N]; // Normalizations are sent upwards

    // assign norm_wire[N] = norm_in;

    // Dummy processing element
    pe_dummy dummy_pe (
        .clk(clk),
        .reset(reset),
        .doProcess(doProcess),
        .data_in(norm_in),
        .data_out(norm_wire[N]) // Connect dummy output to norm_wire start
    );

    // Instantiate PEs
    generate
        for (genvar i = 0; i < N; i++) begin : row
            pe_norm pe (
                .clk(clk),
                .reset(reset),
                .doProcess(doProcess),
                .x_in(data_in[i]),
                .norm_in(norm_wire[i+1]),
                .x_normed_out(data_normed_out[i]),
                .norm_out(norm_wire[i])
            );
        end
    endgenerate

endmodule

`default_nettype wire
