`timescale 1ns / 1ps
`default_nettype none


module pe_norm #(
)(
    input logic clk,
    input logic reset,
    input logic doProcess,
    input real x_in,             // Number to be normalized
    input real norm_in,          // Normalizing factor
    output real x_normed_out,    // Normalized number
    output real norm_out         // Normalizing factor propagated to the next PE
);

    real normed_term;
    always_comb begin
        normed_term = x_in / norm_in;
    end
    
    // Register to store the accumulated partial sum
    always_ff @(posedge clk) begin
        if (reset) begin
            x_normed_out <= '0;
            norm_out <= '0;
        end else if (doProcess) begin
            x_normed_out <= normed_term;
            norm_out <= norm_in; 
        end
    end

endmodule

`default_nettype wire
