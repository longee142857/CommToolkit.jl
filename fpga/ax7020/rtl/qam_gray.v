`timescale 1ns/1ps
`default_nettype none

// 方形 Gray M-QAM 的硬映射和硬判决。M 只接受 16 或 64。
// 这是 Julia QAM(M)，不是 QPSK()。QAM(4) 不在这里：它和 QPSK 几何相同、比特不同。
//
// sym_bits[K-1] 是块里的第一位（MSB）。前 log2(sqrt(M)) 位是 I 路 Gray 标号，
// 其余是 Q 路 Gray 标号。自然 PAM 下标 m 的 Gray 标号是 m xor (m >> 1)，
// 电平是 -(sqrt(M)-1) + 2*m，再乘单位 U。
//   16-QAM：U = round(1024 / sqrt(10)) = 324，电平 ±324、±972。
//   64-QAM：U = round(1024 / sqrt(42)) = 158，电平 ±158、±474、±790、±1106。
// 浮点归一化是除以 sqrt(2*(M-1)/3)，使 Es = 1。这里只取最近整数 U。
//
// 判决门限在相邻电平的正中间。正好落在门限上时归到右侧那一格
// （比较用的是 <，不是 <=）。Julia 的 findmin 在距离完全相等时取星座表里
// 更小的下标，只在这个平局上可能差一格。无噪声的星座点不在门限上。

module qam_mapper #(
    parameter integer M = 16
) (
    input  wire [5:0]         sym_bits,
    output wire signed [15:0] i_sym,
    output wire signed [15:0] q_sym
);

    generate
        if ((M != 16) && (M != 64)) begin : g_bad_m
            invalid_qam_mapper_M_must_be_16_or_64 u_bad();
        end
    endgenerate

    function signed [15:0] map16;
        input [1:0] g;
        begin
            case (g)
                2'b00: map16 = -16'sd972;
                2'b01: map16 = -16'sd324;
                2'b11: map16 =  16'sd324;
                default: map16 = 16'sd972;
            endcase
        end
    endfunction

    function signed [15:0] map64;
        input [2:0] g;
        begin
            case (g)
                3'b000: map64 = -16'sd1106;
                3'b001: map64 = -16'sd790;
                3'b011: map64 = -16'sd474;
                3'b010: map64 = -16'sd158;
                3'b110: map64 =  16'sd158;
                3'b111: map64 =  16'sd474;
                3'b101: map64 =  16'sd790;
                default: map64 = 16'sd1106;
            endcase
        end
    endfunction

    generate
        if (M == 16) begin : g16
            assign i_sym = map16(sym_bits[3:2]);
            assign q_sym = map16(sym_bits[1:0]);
        end else if (M == 64) begin : g64
            assign i_sym = map64(sym_bits[5:3]);
            assign q_sym = map64(sym_bits[2:0]);
        end else begin : g_none
            assign i_sym = 16'sd0;
            assign q_sym = 16'sd0;
        end
    endgenerate

endmodule

module qam_slicer #(
    parameter integer M = 16
) (
    input  wire signed [15:0] i_sym,
    input  wire signed [15:0] q_sym,
    output wire [5:0]         sym_bits
);

    generate
        if ((M != 16) && (M != 64)) begin : g_bad_m
            invalid_qam_slicer_M_must_be_16_or_64 u_bad();
        end
    endgenerate

    function [1:0] slice16;
        input signed [15:0] x;
        begin
            if (x < -16'sd648)
                slice16 = 2'b00;
            else if (x < 16'sd0)
                slice16 = 2'b01;
            else if (x < 16'sd648)
                slice16 = 2'b11;
            else
                slice16 = 2'b10;
        end
    endfunction

    function [2:0] slice64;
        input signed [15:0] x;
        begin
            if (x < -16'sd948)
                slice64 = 3'b000;
            else if (x < -16'sd632)
                slice64 = 3'b001;
            else if (x < -16'sd316)
                slice64 = 3'b011;
            else if (x < 16'sd0)
                slice64 = 3'b010;
            else if (x < 16'sd316)
                slice64 = 3'b110;
            else if (x < 16'sd632)
                slice64 = 3'b111;
            else if (x < 16'sd948)
                slice64 = 3'b101;
            else
                slice64 = 3'b100;
        end
    endfunction

    generate
        if (M == 16) begin : g16
            assign sym_bits = {2'b00, slice16(i_sym), slice16(q_sym)};
        end else if (M == 64) begin : g64
            assign sym_bits = {slice64(i_sym), slice64(q_sym)};
        end else begin : g_none
            assign sym_bits = 6'b000000;
        end
    endgenerate

endmodule

`default_nettype wire
