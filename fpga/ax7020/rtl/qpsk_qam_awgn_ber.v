`timescale 1ns/1ps
`default_nettype none

// 刀2 链路：PRBS → QPSK 或方形 Gray QAM → 两路实 AWGN → 硬判决 → 比特差错。
// MODE = 0 QPSK（每符号 2 比特），1 为 16-QAM（4 比特），2 为 64-QAM（6 比特）。
// 噪声模块与刀1 相同，k 和正交路种子不同。NBITS 必须是每符号比特数的整数倍。
// 一个时钟一个符号。不做法成形、软信息、OFDM、MIMO 或信道码。
module qpsk_qam_awgn_ber #(
    parameter integer MODE  = 0,
    parameter integer NBITS = 16380
) (
    input  wire        clk,
    input  wire        rst,
    input  wire        start,
    input  wire [2:0]  snr_index,
    output wire        busy,
    output wire        done,
    output wire [31:0] n_err,
    output wire [31:0] n_bit,
    output wire [5:0]  tx_bits,
    output wire [5:0]  rx_bits,
    output wire signed [15:0] i_tx,
    output wire signed [15:0] q_tx,
    output wire signed [15:0] i_rx,
    output wire signed [15:0] q_rx
);

    localparam integer K = (MODE == 0) ? 2 : ((MODE == 1) ? 4 : 6);
    localparam integer BPS_OK = (MODE == 0) || (MODE == 1) || (MODE == 2);

    localparam [1:0] S_IDLE = 2'd0;
    localparam [1:0] S_LOAD = 2'd1;
    localparam [1:0] S_RUN  = 2'd2;
    localparam [1:0] S_DONE = 2'd3;

    generate
        if (!BPS_OK) begin : g_mode
            invalid_qpsk_qam_mode u_mode();
        end
        if ((NBITS % K) != 0) begin : g_nbits
            invalid_nbits_not_multiple_of_k u_nbits();
        end
    endgenerate

    reg [1:0]  state;
    reg [1:0]  state_n;
    reg [2:0]  snr_q;
    reg [31:0] n_err_r;
    reg [31:0] n_bit_r;

    wire lfsr_rst = rst || (state == S_IDLE) || (state == S_LOAD);
    wire run_en   = (state == S_RUN);
    wire [7:0] data_byte;

    lfsr_fib #(
        .WIDTH(15),
        .TAP_MASK(15'h4001),
        .STEP(K)
    ) u_data (
        .clk(clk),
        .rst(lfsr_rst),
        .en(run_en),
        .seed(15'h7FFF),
        .bit_out(),
        .sample8(data_byte)
    );

    // data_byte[0] 是这一拍最早的比特，也就是标号的 MSB。
    genvar bi;
    generate
        for (bi = 0; bi < 6; bi = bi + 1) begin : g_pack
            if (bi < K) begin : g_use
                assign tx_bits[bi] = data_byte[(K - 1) - bi];
            end else begin : g_pad
                assign tx_bits[bi] = 1'b0;
            end
        end
    endgenerate

    wire signed [15:0] noise_i;
    wire signed [15:0] noise_q;
    wire [2:0] bps = K;

    awgn_clt #(
        .AMP(1024),
        .SEED_SET(0)
    ) u_ni (
        .clk(clk),
        .rst(lfsr_rst),
        .en(run_en),
        .snr_index(snr_q),
        .bits_per_symbol(bps),
        .noise(noise_i)
    );

    awgn_clt #(
        .AMP(1024),
        .SEED_SET(1)
    ) u_nq (
        .clk(clk),
        .rst(lfsr_rst),
        .en(run_en),
        .snr_index(snr_q),
        .bits_per_symbol(bps),
        .noise(noise_q)
    );

    assign i_rx = i_tx + noise_i;
    assign q_rx = q_tx + noise_q;

    generate
        if (MODE == 0) begin : g_qpsk
            wire [1:0] rx2;
            qpsk_mapper u_map (
                .sym_bits(tx_bits[1:0]),
                .i_sym(i_tx),
                .q_sym(q_tx)
            );
            qpsk_slicer u_slice (
                .i_sym(i_rx),
                .q_sym(q_rx),
                .sym_bits(rx2)
            );
            assign rx_bits = {4'b0000, rx2};
        end else if (MODE == 1) begin : g_qam16
            qam_mapper #(.M(16)) u_map (
                .sym_bits(tx_bits),
                .i_sym(i_tx),
                .q_sym(q_tx)
            );
            qam_slicer #(.M(16)) u_slice (
                .i_sym(i_rx),
                .q_sym(q_rx),
                .sym_bits(rx_bits)
            );
        end else if (MODE == 2) begin : g_qam64
            qam_mapper #(.M(64)) u_map (
                .sym_bits(tx_bits),
                .i_sym(i_tx),
                .q_sym(q_tx)
            );
            qam_slicer #(.M(64)) u_slice (
                .i_sym(i_rx),
                .q_sym(q_rx),
                .sym_bits(rx_bits)
            );
        end else begin : g_tie
            assign i_tx = 16'sd0;
            assign q_tx = 16'sd0;
            assign rx_bits = 6'b000000;
        end
    endgenerate

    // 每项扩成 3 位再相加，避免 1 位加法把一拍里的多个差错折掉。
    wire [2:0] n_err_sym =
        {2'b00, tx_bits[0] ^ rx_bits[0]} +
        {2'b00, tx_bits[1] ^ rx_bits[1]} +
        ((K > 2) ? {2'b00, tx_bits[2] ^ rx_bits[2]} : 3'b000) +
        ((K > 3) ? {2'b00, tx_bits[3] ^ rx_bits[3]} : 3'b000) +
        ((K > 4) ? {2'b00, tx_bits[4] ^ rx_bits[4]} : 3'b000) +
        ((K > 5) ? {2'b00, tx_bits[5] ^ rx_bits[5]} : 3'b000);

    assign busy  = (state == S_LOAD) || (state == S_RUN);
    assign done  = (state == S_DONE);
    assign n_err = n_err_r;
    assign n_bit = n_bit_r;

    always @* begin
        state_n = state;
        case (state)
            S_IDLE: if (start) state_n = S_LOAD;
            S_LOAD: state_n = S_RUN;
            S_RUN:  if (n_bit_r == (NBITS - K)) state_n = S_DONE;
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
                n_bit_r <= n_bit_r + K;
                n_err_r <= n_err_r + n_err_sym;
            end
        end
    end

endmodule

`default_nettype wire
