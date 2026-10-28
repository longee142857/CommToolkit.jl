`timescale 1ns/1ps
`default_nettype none

// 一条可综合的教学链路：PRBS → BPSK → 实 AWGN → 硬判决 → 差错计数。
// 一个符号一个时钟。start 在 IDLE 或 DONE 时被接受；计数满 NBITS 后停在 DONE，
// 寄存器 n_err / n_bit 保持，直到下一次 start。
//
// 数据比特来自 x^15 + x^14 + 1、全 1 种子，与
// lfsr_sequence(prbs_polynomial(15); seed=nothing) 的移位规则相同。
// 噪声不是这条 LFSR 的比特，见 awgn_clt。
module bpsk_awgn_ber #(
    parameter integer NBITS = 16384,
    parameter integer AMP   = 1024
) (
    input  wire               clk,
    input  wire               rst,
    input  wire               start,
    input  wire [2:0]         snr_index,
    output wire               busy,
    output wire               done,
    output wire [31:0]        n_err,
    output wire [31:0]        n_bit,
    output wire               tx_bit,
    output wire               rx_bit,
    output wire signed [15:0] symbol,
    output wire signed [15:0] noise,
    output wire signed [15:0] sample
);

    localparam [1:0] S_IDLE = 2'd0;
    localparam [1:0] S_LOAD = 2'd1;
    localparam [1:0] S_RUN  = 2'd2;
    localparam [1:0] S_DONE = 2'd3;

    reg [1:0]  state;
    reg [1:0]  state_n;
    reg [2:0]  snr_q;
    reg [31:0] n_err_r;
    reg [31:0] n_bit_r;

    wire lfsr_rst = rst || (state == S_IDLE) || (state == S_LOAD);
    wire run_en   = (state == S_RUN);

    lfsr_fib #(
        .WIDTH(15),
        .TAP_MASK(15'h4001), // x^15 + x^14 + 1 的低 15 位：bit 14 | bit 0
        .STEP(1)
    ) u_data (
        .clk(clk),
        .rst(lfsr_rst),
        .en(run_en),
        .seed(15'h7FFF),
        .bit_out(tx_bit),
        .sample8()
    );

    bpsk_mapper #(
        .AMP(AMP)
    ) u_map (
        .bit_in(tx_bit),
        .symbol(symbol)
    );

    awgn_clt #(
        .AMP(AMP)
    ) u_awgn (
        .clk(clk),
        .rst(lfsr_rst),
        .en(run_en),
        .snr_index(snr_q),
        .bits_per_symbol(3'd1),
        .noise(noise)
    );

    assign sample = symbol + noise;

    bpsk_slicer u_slice (
        .sample(sample),
        .bit_out(rx_bit)
    );

    assign busy  = (state == S_LOAD) || (state == S_RUN);
    assign done  = (state == S_DONE);
    assign n_err = n_err_r;
    assign n_bit = n_bit_r;

    always @* begin
        state_n = state;
        case (state)
            S_IDLE: if (start) state_n = S_LOAD;
            S_LOAD: state_n = S_RUN;
            S_RUN:  if (n_bit_r == (NBITS - 1)) state_n = S_DONE;
            S_DONE: if (start) state_n = S_LOAD;
            default: state_n = S_IDLE;
        endcase
    end

    always @(posedge clk) begin
        if (rst) begin
            state   <= S_IDLE;
            snr_q   <= 3'd0;
            n_err_r <= 32'd0;
            n_bit_r <= 32'd0;
        end else begin
            state <= state_n;
            if ((state == S_IDLE || state == S_DONE) && start)
                snr_q <= snr_index;
            if (state == S_LOAD) begin
                n_err_r <= 32'd0;
                n_bit_r <= 32'd0;
            end else if (state == S_RUN) begin
                n_bit_r <= n_bit_r + 32'd1;
                if (tx_bit != rx_bit)
                    n_err_r <= n_err_r + 32'd1;
            end
        end
    end

endmodule

`default_nettype wire
