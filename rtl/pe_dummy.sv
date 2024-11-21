`timescale 1ns / 1ps
`default_nettype none

module pe_dummy #(
)(
    input logic clk,
    input logic reset,
    input logic doProcess,
    input real data_in,             // Input data
    output real data_out         // Output data (same as input)
);
    
    // Register to store the accumulated partial sum
    always_ff @(posedge clk) begin
        if (reset) begin
            data_out <= '0;
        end else if (doProcess) begin
            data_out <= data_in; 
        end
    end

endmodule

`default_nettype wire
