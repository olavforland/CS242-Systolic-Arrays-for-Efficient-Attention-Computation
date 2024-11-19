`timescale 1ns / 1ps
`default_nettype none

module pe #(
)(
    input logic clk,
    input logic reset,
    input logic doProcess,
    input logic weight_load_enable,
    input real weight_in,
    input real data_in,          // Input data from the bottom
    input real partial_sum_in,    // Partial sum from the left
    output real data_out,        // Output data to the top
    output real partial_sum_out   // Output partial sum 
);
    // Internal weight register
    real weight;

    // Load weight on weight_load_enable signal
    always_ff @(posedge clk) begin
        if (reset) begin
            weight <= '0;
        end else if (weight_load_enable) begin
            weight <= weight_in;
        end
    end

    // Register data input to data output
    always_ff @(posedge clk) begin
        if (reset) begin
            data_out <= '0;
        end else if (doProcess) begin
            data_out <= data_in;
        end
    end

    // Multiply-Accumulate operation
    real mul_result;
    real add_result;

    always_comb begin
        mul_result = data_in * weight;
        // Extend mul_result to ACC_WIDTH if necessary
        add_result = mul_result + partial_sum_in;
    end

    // Register to store the accumulated partial sum
    always_ff @(posedge clk) begin
        if (reset) begin
            partial_sum_out <= '0;
        end else if (doProcess) begin
            partial_sum_out <= add_result; 
        end
    end

endmodule

`default_nettype wire
