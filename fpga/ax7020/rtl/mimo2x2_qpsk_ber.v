`timescale 1ns/1ps
`default_nettype none

// 固定 2×2 信道上的 QPSK。理想 CSI。MODE=0 是 ZF（σ²=0），MODE=1 是 MMSE。
// H 的整数刻度是 1024 = 1。默认矩阵 [[1024, 256], [128, 1024]]，实数，列不满正交，
// 所以两路之间有干扰。不是 `mimo_channel` 抽出来的瑞利矩阵，也不对 H 做平均。
//
// 每个接收天线上的实噪声标准差是刀2 QPSK 刻度 σ 再乘 1024，和「1.0」的刻度对齐。
// MMSE 加到 Gram 上的是 2σ²，不是 y 上那个已经乘过 1024 的样本方差。
// 差一个 1024² 是为了让 64 位行列式装得下；因此这条 MMSE 比把真实 y 方差代进去更靠近 ZF。
// 见 README。snr_index >= 6 时 σ=0，两种模式都变成无噪声 ZF。
module mimo2x2_qpsk_ber #(
    parameter integer MODE  = 0,
    parameter integer NVEC  = 256,
    parameter signed [15:0] H00_RE = 16'sd1024,
    parameter signed [15:0] H00_IM = 16'sd0,
    parameter signed [15:0] H01_RE = 16'sd256,
    parameter signed [15:0] H01_IM = 16'sd0,
    parameter signed [15:0] H10_RE = 16'sd128,
    parameter signed [15:0] H10_IM = 16'sd0,
    parameter signed [15:0] H11_RE = 16'sd1024,
    parameter signed [15:0] H11_IM = 16'sd0
) (
    input  wire        clk,
    input  wire        rst,
    input  wire        start,
    input  wire [2:0]  snr_index,
    output wire        busy,
    output wire        done,
    output wire [31:0] n_err,
    output wire [31:0] n_bit
);

    localparam [2:0] S_IDLE  = 3'd0;
    localparam [2:0] S_LOAD  = 3'd1;
    localparam [2:0] S_SYMB  = 3'd2;
    localparam [2:0] S_N0    = 3'd3;
    localparam [2:0] S_N1    = 3'd4;
    localparam [2:0] S_COUNT = 3'd5;
    localparam [2:0] S_DONE  = 3'd6;

    reg [2:0]  state;
    reg [2:0]  snr_q;
    reg [15:0] vec;
    reg [31:0] n_err_r;
    reg [31:0] n_bit_r;
    reg [3:0]  bits_q;
    reg signed [31:0] n0_re, n0_im, n1_re, n1_im;

    wire lfsr_rst = rst || (state == S_IDLE) || (state == S_LOAD);
    wire data_en  = (state == S_SYMB);
    wire noise_en = (state == S_N0) || (state == S_N1);

    wire [7:0] data_byte;
    lfsr_fib #(
        .WIDTH(15),
        .TAP_MASK(15'h4001),
        .STEP(4)
    ) u_data (
        .clk(clk),
        .rst(lfsr_rst),
        .en(data_en),
        .seed(15'h7FFF),
        .bit_out(),
        .sample8(data_byte)
    );

    // [0] 最早，是第 0 路的 b0。
    wire [1:0] bits0 = {data_byte[0], data_byte[1]};
    wire [1:0] bits1 = {data_byte[2], data_byte[3]};
    wire signed [15:0] x0_i, x0_q, x1_i, x1_q;
    qpsk_mapper u_x0 (.sym_bits(bits0), .i_sym(x0_i), .q_sym(x0_q));
    qpsk_mapper u_x1 (.sym_bits(bits1), .i_sym(x1_i), .q_sym(x1_q));

    reg signed [15:0] sx0_i, sx0_q, sx1_i, sx1_q;

    wire signed [15:0] noise_i, noise_q;
    awgn_clt #(.AMP(1024), .SEED_SET(0)) u_ni (
        .clk(clk), .rst(lfsr_rst), .en(noise_en),
        .snr_index(snr_q), .bits_per_symbol(3'd2), .noise(noise_i)
    );
    awgn_clt #(.AMP(1024), .SEED_SET(1)) u_nq (
        .clk(clk), .rst(lfsr_rst), .en(noise_en),
        .snr_index(snr_q), .bits_per_symbol(3'd2), .noise(noise_q)
    );

    function signed [31:0] ymul_re;
        input signed [15:0] hr, hi, xr, xi;
        begin
            ymul_re = hr * xr - hi * xi;
        end
    endfunction
    function signed [31:0] ymul_im;
        input signed [15:0] hr, hi, xr, xi;
        begin
            ymul_im = hr * xi + hi * xr;
        end
    endfunction

    wire signed [31:0] y0_re =
        ymul_re(H00_RE, H00_IM, sx0_i, sx0_q) +
        ymul_re(H01_RE, H01_IM, sx1_i, sx1_q) + n0_re;
    wire signed [31:0] y0_im =
        ymul_im(H00_RE, H00_IM, sx0_i, sx0_q) +
        ymul_im(H01_RE, H01_IM, sx1_i, sx1_q) + n0_im;
    wire signed [31:0] y1_re =
        ymul_re(H10_RE, H10_IM, sx0_i, sx0_q) +
        ymul_re(H11_RE, H11_IM, sx1_i, sx1_q) + n1_re;
    wire signed [31:0] y1_im =
        ymul_im(H10_RE, H10_IM, sx0_i, sx0_q) +
        ymul_im(H11_RE, H11_IM, sx1_i, sx1_q) + n1_im;

    function [31:0] nvar_of;
        input [2:0] idx;
        input       mmse;
        reg [15:0] s;
        begin
            if (!mmse) begin
                nvar_of = 32'd0;
            end else begin
                case (idx)
                    3'd0: s = 16'd512;
                    3'd1: s = 16'd407;
                    3'd2: s = 16'd323;
                    3'd3: s = 16'd257;
                    3'd4: s = 16'd204;
                    3'd5: s = 16'd162;
                    default: s = 16'd0;
                endcase
                nvar_of = (s * s) << 1;
            end
        end
    endfunction

    wire [31:0] noise_var = nvar_of(snr_q, MODE != 0);
    wire signed [15:0] xh0_i, xh0_q, xh1_i, xh1_q;
    wire singular;

    mimo2x2_detect u_det (
        .h00_re(H00_RE), .h00_im(H00_IM),
        .h01_re(H01_RE), .h01_im(H01_IM),
        .h10_re(H10_RE), .h10_im(H10_IM),
        .h11_re(H11_RE), .h11_im(H11_IM),
        .y0_re(y0_re), .y0_im(y0_im),
        .y1_re(y1_re), .y1_im(y1_im),
        .noise_var(noise_var),
        .x0_re(xh0_i), .x0_im(xh0_q),
        .x1_re(xh1_i), .x1_im(xh1_q),
        .singular(singular)
    );

    wire [1:0] hat0, hat1;
    qpsk_slicer u_s0 (.i_sym(xh0_i), .q_sym(xh0_q), .sym_bits(hat0));
    qpsk_slicer u_s1 (.i_sym(xh1_i), .q_sym(xh1_q), .sym_bits(hat1));

    // 噪声从 σ 放大 1024 倍，对齐 H 的「1.0」。
    wire signed [31:0] nbig_i = {{6{noise_i[15]}}, noise_i, 10'b0};
    wire signed [31:0] nbig_q = {{6{noise_q[15]}}, noise_q, 10'b0};

    assign busy  = (state != S_IDLE) && (state != S_DONE);
    assign done  = (state == S_DONE);
    assign n_err = n_err_r;
    assign n_bit = n_bit_r;

    always @(posedge clk) begin
        if (rst) begin
            state   <= S_IDLE;
            snr_q   <= 3'd0;
            vec     <= 16'd0;
            n_err_r <= 32'd0;
            n_bit_r <= 32'd0;
            bits_q  <= 4'd0;
            sx0_i <= 0; sx0_q <= 0; sx1_i <= 0; sx1_q <= 0;
            n0_re <= 0; n0_im <= 0; n1_re <= 0; n1_im <= 0;
        end else begin
            case (state)
                S_IDLE: begin
                    if (start) begin
                        state   <= S_LOAD;
                        snr_q   <= snr_index;
                        vec     <= 16'd0;
                        n_err_r <= 32'd0;
                        n_bit_r <= 32'd0;
                    end
                end
                S_LOAD: state <= S_SYMB;
                S_SYMB: begin
                    bits_q <= {data_byte[0], data_byte[1], data_byte[2], data_byte[3]};
                    sx0_i  <= x0_i;
                    sx0_q  <= x0_q;
                    sx1_i  <= x1_i;
                    sx1_q  <= x1_q;
                    state  <= S_N0;
                end
                S_N0: begin
                    n0_re <= nbig_i;
                    n0_im <= nbig_q;
                    state <= S_N1;
                end
                S_N1: begin
                    n1_re <= nbig_i;
                    n1_im <= nbig_q;
                    state <= S_COUNT;
                end
                S_COUNT: begin
                    n_bit_r <= n_bit_r + 32'd4;
                    n_err_r <= n_err_r
                        + {31'd0, singular | (hat0[1] != bits_q[3])}
                        + {31'd0, singular | (hat0[0] != bits_q[2])}
                        + {31'd0, singular | (hat1[1] != bits_q[1])}
                        + {31'd0, singular | (hat1[0] != bits_q[0])};
                    if (vec == (NVEC - 1))
                        state <= S_DONE;
                    else begin
                        vec   <= vec + 16'd1;
                        state <= S_SYMB;
                    end
                end
                S_DONE: begin
                    if (start) begin
                        state   <= S_LOAD;
                        snr_q   <= snr_index;
                        vec     <= 16'd0;
                        n_err_r <= 32'd0;
                        n_bit_r <= 32'd0;
                    end
                end
                default: state <= S_IDLE;
            endcase
        end
    end

endmodule

`default_nettype wire
