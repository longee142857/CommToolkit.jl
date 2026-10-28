`timescale 1ns/1ps
`default_nettype none

// 信息比特 → (7, 5) 编码（含 2 个 0 尾比特）→ BPSK → 两路实 AWGN → 硬判决
// → 硬 Viterbi → 信息比特差错。一个编码步一个时钟，两路噪声对应 (g1, g2)。
//
// snr_index 是信息比特 Eb/N0。名义码率 1/2、编码后 BPSK 的 Es = 1，所以
// Eb = 2，和 Julia 测试里 awgn(..., energy_per_bit=2) 一样。尾比特使块码率
// 变成 n/(2n+4)，这个损失没有再折进噪声。SCALE_MODE=1 选的就是这张刻度。
// 无噪声用 snr_index >= 6。差错只数信息比特，不含尾比特。
module conv_bpsk_awgn_ber #(
    parameter integer NINFO   = 64,
    parameter integer NBLOCKS = 32
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

    localparam integer NSYM = NINFO + 2;

    localparam [2:0] S_IDLE = 3'd0;
    localparam [2:0] S_LOAD = 3'd1;
    localparam [2:0] S_BLK  = 3'd2;
    localparam [2:0] S_ENC  = 3'd3;
    localparam [2:0] S_WAIT = 3'd4;
    localparam [2:0] S_DONE = 3'd5;

    reg [2:0]  state;
    reg [2:0]  snr_q;
    reg [9:0]  sym;
    reg [15:0] blk;
    reg [31:0] n_err_r;
    reg [31:0] n_bit_r;
    reg        info_mem [0:NINFO-1];

    wire lfsr_rst = rst || (state == S_IDLE) || (state == S_LOAD);
    wire data_en  = (state == S_ENC) && (sym < NINFO);
    wire noise_en = (state == S_ENC);
    wire enc_rst  = rst || (state == S_BLK);
    wire enc_en   = (state == S_ENC);
    wire tail     = (sym >= NINFO);

    wire data_bit;
    lfsr_fib #(
        .WIDTH(15),
        .TAP_MASK(15'h4001),
        .STEP(1)
    ) u_data (
        .clk(clk),
        .rst(lfsr_rst),
        .en(data_en),
        .seed(15'h7FFF),
        .bit_out(data_bit),
        .sample8()
    );

    wire       bit_in = tail ? 1'b0 : data_bit;
    wire [1:0] coded;
    wire [1:0] enc_state;
    conv75_encode u_enc (
        .clk(clk),
        .rst(enc_rst),
        .en(enc_en),
        .bit_in(bit_in),
        .coded(coded),
        .state(enc_state)
    );

    wire signed [15:0] sym_g1, sym_g2, noise_g1, noise_g2;
    bpsk_mapper u_m1 (
        .bit_in(coded[1]),
        .symbol(sym_g1)
    );
    bpsk_mapper u_m0 (
        .bit_in(coded[0]),
        .symbol(sym_g2)
    );

    awgn_clt #(
        .AMP(1024),
        .SEED_SET(0),
        .SCALE_MODE(1)
    ) u_n1 (
        .clk(clk),
        .rst(lfsr_rst),
        .en(noise_en),
        .snr_index(snr_q),
        .bits_per_symbol(3'd1),
        .noise(noise_g1)
    );
    awgn_clt #(
        .AMP(1024),
        .SEED_SET(1),
        .SCALE_MODE(1)
    ) u_n0 (
        .clk(clk),
        .rst(lfsr_rst),
        .en(noise_en),
        .snr_index(snr_q),
        .bits_per_symbol(3'd1),
        .noise(noise_g2)
    );

    wire hard_g1, hard_g0;
    bpsk_slicer u_s1 (
        .sample(sym_g1 + noise_g1),
        .bit_out(hard_g1)
    );
    bpsk_slicer u_s0 (
        .sample(sym_g2 + noise_g2),
        .bit_out(hard_g0)
    );

    wire        dec_start = (state == S_BLK);
    wire        rx_valid  = (state == S_ENC);
    wire        dec_busy, dec_done, info_valid, info_bit;
    wire [9:0]  info_idx;
    viterbi_hard75 #(
        .NINFO(NINFO)
    ) u_dec (
        .clk(clk),
        .rst(rst),
        .start(dec_start),
        .rx_valid(rx_valid),
        .rx_pair({hard_g1, hard_g0}),
        .busy(dec_busy),
        .done(dec_done),
        .info_valid(info_valid),
        .info_bit(info_bit),
        .info_idx(info_idx)
    );

    assign busy  = (state != S_IDLE) && (state != S_DONE);
    assign done  = (state == S_DONE);
    assign n_err = n_err_r;
    assign n_bit = n_bit_r;

    always @(posedge clk) begin
        if (rst) begin
            state   <= S_IDLE;
            snr_q   <= 3'd0;
            sym     <= 10'd0;
            blk     <= 16'd0;
            n_err_r <= 32'd0;
            n_bit_r <= 32'd0;
        end else begin
            if (info_valid) begin
                n_bit_r <= n_bit_r + 32'd1;
                if (info_bit != info_mem[info_idx])
                    n_err_r <= n_err_r + 32'd1;
            end
            case (state)
                S_IDLE: begin
                    if (start) begin
                        state   <= S_LOAD;
                        snr_q   <= snr_index;
                        blk     <= 16'd0;
                        n_err_r <= 32'd0;
                        n_bit_r <= 32'd0;
                    end
                end
                S_LOAD: begin
                    state <= S_BLK;
                    sym   <= 10'd0;
                end
                S_BLK: begin
                    state <= S_ENC;
                    sym   <= 10'd0;
                end
                S_ENC: begin
                    if (sym < NINFO)
                        info_mem[sym] <= data_bit;
                    if (sym == (NSYM - 1))
                        state <= S_WAIT;
                    else
                        sym <= sym + 10'd1;
                end
                S_WAIT: begin
                    if (dec_done) begin
                        if (blk == (NBLOCKS - 1))
                            state <= S_DONE;
                        else begin
                            blk   <= blk + 16'd1;
                            sym   <= 10'd0;
                            state <= S_BLK;
                        end
                    end
                end
                S_DONE: begin
                    if (start) begin
                        state   <= S_LOAD;
                        snr_q   <= snr_index;
                        blk     <= 16'd0;
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
