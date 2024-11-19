`timescale 1ns / 1ps
`default_nettype none

module pe_exp #(
)(
    input logic clk,
    input logic reset,
    input logic doProcess,
    input real k_factorial,               // Factorial term of current iteration
    input real x_in,                // Number to be exponentiated
    input real x_pow_in,            // Number to be exponentiated to the power of k
    input real taylor_approx_in,    // Cumulative sum of the Taylor series
    output real x_out,               // Number to be exponentiated
    output real x_pow_out,           // Number to be exponentiated to the power of k+1
    output real taylor_approx_out    // Cumulative sum after adding the Taylor series term
);

    // Multiply-Accumulate operation
    real exp_term, taylor_term, taylor_sum;
    
    always_comb begin
        // Increment the power of x
        exp_term = x_pow_in * x_in;
        // Calculate the factorial term
        taylor_term = exp_term / k_factorial;
        // Add the term to the cumulative sum
        taylor_sum = taylor_approx_in + taylor_term;
    end

    // Register to store the accumulated partial sum
    always_ff @(posedge clk) begin
        if (reset) begin
            x_out <= '0;
            x_pow_out <= '0;
            taylor_approx_out <= '0;
        end else if (doProcess) begin
            x_out <= x_in;
            x_pow_out <= exp_term;
            taylor_approx_out <= taylor_sum; 
        end
    end
endmodule

`default_nettype wire
