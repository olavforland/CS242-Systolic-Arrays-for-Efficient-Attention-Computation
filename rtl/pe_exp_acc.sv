`timescale 1ns / 1ps
`default_nettype none


module pe_exp_acc #(
)(
    input logic clk,
    input logic reset,
    input logic doProcess,
    input real exp_in,                // Final Taylor series approximation
    input real exp_sum_in,            // Cumulative sum of the Taylor series terms
    output real exp_out,               // Final Taylor series approximation (same as exp_in)
    output real exp_sum_out            // exp_in + exp_sum_in
);
    real cumulative_sum;
    always begin
        cumulative_sum = exp_in + exp_sum_in;
    end
    
    // Register to store the accumulated partial sum
    always_ff @(posedge clk) begin
        if (reset) begin
            exp_out <= '0;
            exp_sum_out <= '0;
        end else if (doProcess) begin
            exp_out <= exp_in;
            exp_sum_out <= cumulative_sum; 
        end
    end

endmodule

`default_nettype wire
