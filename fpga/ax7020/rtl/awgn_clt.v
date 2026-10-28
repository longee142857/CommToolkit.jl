`timescale 1ns/1ps
`default_nettype none

// 教学级实 AWGN 近似。不是 Box–Muller，也不是高斯 ROM。
//
// 12 路独立种子的 31 级 Fibonacci LFSR（多项式 x^31 + x^28 + 1，与
// prbs_polynomial(31) 的连接多项式相同）各取当前状态低 8 位，当作
// {0,...,255} 上的离散均匀数，相加后减去均值 12*127.5 = 1530。
// 和的标准差是 sqrt(12 * (256^2-1)/12) = sqrt(65535) ≈ 255.998。
// 再乘 Eb/N0 刻度并算术右移 8 位，得到近似 N(0, σ^2) 的整数噪声：
//
//   σ = AMP / sqrt(2 * 10^(Eb/N0_dB / 10))
//   MULT = round(σ)，SHIFT = 8
//   noise = floor(centered * MULT / 256)   // 负数朝 -∞ 截断
//
// snr_index: 0..5 → 0,2,4,6,8,10 dB；6 和 7 → MULT = 0（无噪声）。
// 分布是离散均匀和（Irwin–Hall 型），两边被截断，不是 Julia 的 randn。
// 因此误码率只要求随 Eb/N0 下降，不与浮点仿真或 Q 函数做数值契约。
module awgn_clt #(
    parameter integer AMP = 1024
) (
    input  wire               clk,
    input  wire               rst,
    input  wire               en,
    input  wire [2:0]         snr_index,
    output wire signed [15:0] noise
);

    localparam integer TAP31 = 31'h1000_0001; // bit 28 | bit 0

    // 上面的 MULT 表只按 AMP = 1024 算过。换幅度必须重算整张表，否则 Eb/N0 是错的。
    generate
        if (AMP != 1024) begin : g_amp
            invalid_awgn_amp_only_1024_is_scaled u_scale_not_retuned();
        end
    endgenerate

    function [15:0] scale_of;
        input [2:0] idx;
        begin
            // round(AMP / sqrt(2 * 10^(dB/10)))，AMP = 1024。
            // 0 dB 724, 2 dB 575, 4 dB 457, 6 dB 363, 8 dB 288, 10 dB 229。
            case (idx)
                3'd0: scale_of = 16'd724;
                3'd1: scale_of = 16'd575;
                3'd2: scale_of = 16'd457;
                3'd3: scale_of = 16'd363;
                3'd4: scale_of = 16'd288;
                3'd5: scale_of = 16'd229;
                default: scale_of = 16'd0;
            endcase
        end
    endfunction

    // 12 个非零 31 位种子。只决定可复现序列，不是 Julia RNG 的状态。
    function [30:0] seed_of;
        input integer idx;
        begin
            case (idx)
                0:  seed_of = 31'h2E4B17C3;
                1:  seed_of = 31'h51A60D29;
                2:  seed_of = 31'h0C0DE116;
                3:  seed_of = 31'h6A5A1235;
                4:  seed_of = 31'h13579BDF;
                5:  seed_of = 31'h2468ACE0;
                6:  seed_of = 31'h0F1E2D3C;
                7:  seed_of = 31'h7FFFFFFE;
                8:  seed_of = 31'h11111111;
                9:  seed_of = 31'h40000001;
                10: seed_of = 31'h00ABCDEF;
                default: seed_of = 31'h55AA33CC;
            endcase
        end
    endfunction

    wire [7:0] ubyte [0:11];
    wire       unused_bit [0:11];

    genvar gi;
    generate
        for (gi = 0; gi < 12; gi = gi + 1) begin : g_noise
            lfsr_fib #(
                .WIDTH(31),
                .TAP_MASK(TAP31),
                .STEP(8)
            ) u_lfsr (
                .clk(clk),
                .rst(rst),
                .en(en),
                .seed(seed_of(gi)),
                .bit_out(unused_bit[gi]),
                .sample8(ubyte[gi])
            );
        end
    endgenerate

    // 每项先扩到 12 位。12×255 = 3060，加总不会溢出。
    wire [11:0] sum_u =
        {4'd0, ubyte[0]}  + {4'd0, ubyte[1]}  + {4'd0, ubyte[2]}  +
        {4'd0, ubyte[3]}  + {4'd0, ubyte[4]}  + {4'd0, ubyte[5]}  +
        {4'd0, ubyte[6]}  + {4'd0, ubyte[7]}  + {4'd0, ubyte[8]}  +
        {4'd0, ubyte[9]}  + {4'd0, ubyte[10]} + {4'd0, ubyte[11]};

    wire signed [12:0] centered = $signed({1'b0, sum_u}) - 13'sd1530;
    wire [15:0]        mult     = scale_of(snr_index);
    wire signed [16:0] mult_s   = $signed({1'b0, mult});
    wire signed [29:0] prod     = centered * mult_s;
    wire signed [29:0] shifted  = prod >>> 8;

    assign noise = shifted[15:0];

endmodule

`default_nettype wire
