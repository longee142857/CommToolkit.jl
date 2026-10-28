`default_nettype none

// BPSK 映射，概念上对应 modulate(bits, BPSK()) = 1 - 2*bits。
// 0 → +AMP，1 → -AMP。AMP 是整数幅度，不是 Julia 的浮点 +1/-1。
// AMP 必须为正，且 |AMP| 与最大噪声之和要落在有符号 16 位里（见 README）。
module bpsk_mapper #(
    parameter integer AMP = 1024
) (
    input  wire              bit_in,
    output wire signed [15:0] symbol
);

    localparam signed [15:0] AMP_Q = AMP;

    assign symbol = bit_in ? -AMP_Q : AMP_Q;

endmodule

`default_nettype wire
