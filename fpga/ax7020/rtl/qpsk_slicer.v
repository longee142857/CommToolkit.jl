`timescale 1ns/1ps
`default_nettype none

// QPSK 硬判决，对应 demodulate(::QPSK)：
// b0 = (I < 0)，b1 = (Q < 0)。0 落在正半边，判成比特 0。
module qpsk_slicer (
    input  wire signed [15:0] i_sym,
    input  wire signed [15:0] q_sym,
    output wire [1:0]         sym_bits
);

    assign sym_bits[1] = (i_sym < 16'sd0);
    assign sym_bits[0] = (q_sym < 16'sd0);

endmodule

`default_nettype wire
