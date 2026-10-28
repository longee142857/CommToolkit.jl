`timescale 1ns/1ps
`default_nettype none

// Julia QPSK() 的硬映射，不是 QAM(4)，也不是 PSK(4)。
// sym_bits[1] 是先送入的 b0，sym_bits[0] 是 b1。
//   00 → (+A, +A)    01 → (+A, -A)
//   10 → (-A, +A)    11 → (-A, -A)
// A = round(1024 / sqrt(2)) = 724。浮点星座是 (±1 ± j) / sqrt(2)，
// 这里用整数近似单位符号能量，不声称与 modulate 的浮点坐标 bit-exact。
module qpsk_mapper (
    input  wire [1:0]         sym_bits,
    output wire signed [15:0] i_sym,
    output wire signed [15:0] q_sym
);

    localparam signed [15:0] A = 16'sd724;

    assign i_sym = sym_bits[1] ? -A : A;
    assign q_sym = sym_bits[0] ? -A : A;

endmodule

`default_nettype wire
