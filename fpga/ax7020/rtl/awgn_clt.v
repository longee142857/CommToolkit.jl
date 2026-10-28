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
//   σ = AMP / sqrt(2 * k * 10^(Eb/N0_dB / 10))
//   MULT = round(σ)，SHIFT = 8
//   noise = floor(centered * MULT / 256)   // 负数朝 -∞ 截断
//
// k = bits_per_symbol。k = 1 就是刀1 的 BPSK 表（Es = Eb）。
// k = 2、4、6 给 QPSK、16-QAM、64-QAM，对应 Julia
// awgn(syms, snr_db, k) 在 Es = 1 时每个实维度的方差 1/(2*k*ρ)。
// snr_index: 0..5 → 0,2,4,6,8,10 dB；其它 → MULT = 0（无噪声）。
// 其它 k 也得到 MULT = 0，调用方不要靠这个来“关闭噪声”以外的用途。
//
// SEED_SET = 0 是刀1 的 12 个种子。SEED_SET = 1 是另一组，给正交路用，
// 这样 I/Q 两路噪声不会是同一条序列。分布仍是截断的均匀和，不是 randn。
module awgn_clt #(
    parameter integer AMP      = 1024,
    parameter integer SEED_SET = 0
) (
    input  wire               clk,
    input  wire               rst,
    input  wire               en,
    input  wire [2:0]         snr_index,
    input  wire [2:0]         bits_per_symbol,
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
        input [2:0] bps;
        begin
            // round(1024 / sqrt(2 * k * 10^(dB/10)))。k = 1 的六个数与刀1 相同。
            if (bps == 3'd1) begin
                case (idx)
                    3'd0: scale_of = 16'd724;
                    3'd1: scale_of = 16'd575;
                    3'd2: scale_of = 16'd457;
                    3'd3: scale_of = 16'd363;
                    3'd4: scale_of = 16'd288;
                    3'd5: scale_of = 16'd229;
                    default: scale_of = 16'd0;
                endcase
            end else if (bps == 3'd2) begin
                case (idx)
                    3'd0: scale_of = 16'd512;
                    3'd1: scale_of = 16'd407;
                    3'd2: scale_of = 16'd323;
                    3'd3: scale_of = 16'd257;
                    3'd4: scale_of = 16'd204;
                    3'd5: scale_of = 16'd162;
                    default: scale_of = 16'd0;
                endcase
            end else if (bps == 3'd4) begin
                case (idx)
                    3'd0: scale_of = 16'd362;
                    3'd1: scale_of = 16'd288;
                    3'd2: scale_of = 16'd228;
                    3'd3: scale_of = 16'd181;
                    3'd4: scale_of = 16'd144;
                    3'd5: scale_of = 16'd114;
                    default: scale_of = 16'd0;
                endcase
            end else if (bps == 3'd6) begin
                case (idx)
                    3'd0: scale_of = 16'd296;
                    3'd1: scale_of = 16'd235;
                    3'd2: scale_of = 16'd187;
                    3'd3: scale_of = 16'd148;
                    3'd4: scale_of = 16'd118;
                    3'd5: scale_of = 16'd93;
                    default: scale_of = 16'd0;
                endcase
            end else begin
                scale_of = 16'd0;
            end
        end
    endfunction

    // 12 个非零 31 位种子。只决定可复现序列，不是 Julia RNG 的状态。
    function [30:0] seed_of;
        input integer idx;
        begin
            if (SEED_SET == 0) begin
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
            end else begin
                case (idx)
                    0:  seed_of = 31'h1A2B3C4D;
                    1:  seed_of = 31'h5F0E1D2C;
                    2:  seed_of = 31'h3B4C5D6E;
                    3:  seed_of = 31'h01234567;
                    4:  seed_of = 31'h6E5D4C3B;
                    5:  seed_of = 31'h2BADF00D;
                    6:  seed_of = 31'h48ACE135;
                    7:  seed_of = 31'h0000F00D;
                    8:  seed_of = 31'h35C4A1F0;
                    9:  seed_of = 31'h19D28B63;
                    10: seed_of = 31'h3F000001;
                    default: seed_of = 31'h6D5C4B3A;
                endcase
            end
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
    wire [15:0]        mult     = scale_of(snr_index, bits_per_symbol);
    wire signed [16:0] mult_s   = $signed({1'b0, mult});
    wire signed [29:0] prod     = centered * mult_s;
    wire signed [29:0] shifted  = prod >>> 8;

    assign noise = shifted[15:0];

endmodule

`default_nettype wire
