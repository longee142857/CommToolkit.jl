`timescale 1ns/1ps
`default_nettype none

// 教学 MCS / Chase 闭环。不是 NR 调度器，不是增量冗余，没有 BICM。
//
// 三档，snr_index 的 0,2,4,6,8 dB 是 Es/N0（噪声用刀1 的 k=1 表，energy_per_bit=1）：
//   0  未编码 BPSK     16 信息比特 + CRC-16，32 个符号。snr_index 0、1 自动选它。
//   1  (7,5) 母码+BPSK  同上 32 比特再加 2 个尾比特，68 个符号。snr_index 2、3。
//   2  未编码 QPSK      32 比特，16 个符号。snr_index ≥ 4，含无噪声的 6、7。
// mcs_force = 3 才自动选；0、1、2 强制该档。门限是这三档自己的，不是
// teaching_amc_thresholds()。
//
// 重传同一份符号。硬 Chase：每次硬判决记 +1（比特 0）或 -1（比特 1），求和后
// 负数判 1，0 和正数判 1 的反面（判 0）。不是 LLR 相加。CRC 过了就停，否则直到 max_tx。
// 卷积档先把合并后的硬比特送进硬 Viterbi，再查 CRC。
//
// 顶层是可仿真、也按可综合 Verilog 写的状态机，没有跑过 Vivado。
module link_mcs_harq #(
    parameter integer NBLOCKS = 8
) (
    input  wire        clk,
    input  wire        rst,
    input  wire        start,
    input  wire [2:0]  snr_index,
    input  wire [1:0]  mcs_force,
    input  wire [2:0]  max_tx,
    output wire        busy,
    output wire        done,
    output wire [1:0]  mcs_used,
    output wire [15:0] n_ack,
    output wire [15:0] n_fail,
    output wire [31:0] n_uses,
    output wire [31:0] n_good
);

    // CRC 余数在「送入最后一比特」的下一拍才稳定，所以 CRCIN / CCHK
    // 和读取余数的状态之间各空一拍（CRCDONE、CWAIT）。
    localparam [4:0] S_IDLE    = 5'd0;
    localparam [4:0] S_LOAD    = 5'd1;
    localparam [4:0] S_PAY     = 5'd2;
    localparam [4:0] S_CRC0    = 5'd3;
    localparam [4:0] S_CRCIN   = 5'd4;
    localparam [4:0] S_CRCDONE = 5'd5;
    localparam [4:0] S_APPEND  = 5'd6;
    localparam [4:0] S_ENC0    = 5'd7;
    localparam [4:0] S_ENC     = 5'd8;
    localparam [4:0] S_TX      = 5'd9;
    localparam [4:0] S_VARM    = 5'd10;
    localparam [4:0] S_VFEED   = 5'd11;
    localparam [4:0] S_VWAIT   = 5'd12;
    localparam [4:0] S_C0      = 5'd13;
    localparam [4:0] S_CCHK    = 5'd14;
    localparam [4:0] S_CWAIT   = 5'd15;
    localparam [4:0] S_NEXT    = 5'd16;
    localparam [4:0] S_DONE    = 5'd17;

    reg [4:0]  state;
    reg [2:0]  snr_q;
    reg [1:0]  mcs_q;
    reg [2:0]  max_tx_q;
    reg [15:0] blk;
    reg [6:0]  idx;
    reg [2:0]  tx_try;
    reg        saw_low;
    reg [15:0] ack_r, fail_r;
    reg [31:0] uses_r, good_r;

    reg        pay   [0:15];
    reg        word  [0:31];
    reg        chan  [0:67];
    reg        decb  [0:31];
    reg signed [3:0] chase [0:67];

    wire lfsr_rst = rst || (state == S_IDLE) || (state == S_LOAD);
    wire data_en  = (state == S_PAY);
    wire noise_en = (state == S_TX);

    wire bit_lfsr;
    lfsr_fib #(
        .WIDTH(15), .TAP_MASK(15'h4001), .STEP(1)
    ) u_pay (
        .clk(clk), .rst(lfsr_rst), .en(data_en), .seed(15'h7FFF),
        .bit_out(bit_lfsr), .sample8()
    );

    reg        crc_rst, crc_en, crc_bit;
    wire [15:0] crc_reg;
    crc16_ccitt u_crc (
        .clk(clk), .rst(crc_rst), .en(crc_en), .bit_in(crc_bit), .crc(crc_reg)
    );

    reg        enc_rst, enc_en, enc_bit;
    wire [1:0] enc_pair;
    conv75_encode u_enc (
        .clk(clk), .rst(enc_rst), .en(enc_en), .bit_in(enc_bit),
        .coded(enc_pair), .state()
    );

    wire signed [15:0] noise_i, noise_q;
    awgn_clt #(.AMP(1024), .SEED_SET(0)) u_ni (
        .clk(clk), .rst(lfsr_rst), .en(noise_en),
        .snr_index(snr_q), .bits_per_symbol(3'd1), .noise(noise_i)
    );
    awgn_clt #(.AMP(1024), .SEED_SET(1)) u_nq (
        .clk(clk), .rst(lfsr_rst), .en(noise_en),
        .snr_index(snr_q), .bits_per_symbol(3'd1), .noise(noise_q)
    );

    wire [6:0] nsym = (mcs_q == 2'd1) ? 7'd68 : (mcs_q == 2'd2) ? 7'd16 : 7'd32;

    wire        bpsk_bit = (mcs_q == 2'd1) ? chan[idx] : word[idx];
    wire signed [15:0] bpsk_sym;
    bpsk_mapper u_bmap (.bit_in(bpsk_bit), .symbol(bpsk_sym));
    wire bpsk_hat;
    bpsk_slicer u_bsl (.sample(bpsk_sym + noise_i), .bit_out(bpsk_hat));

    wire [1:0] qbits = {word[{idx[3:0], 1'b0}], word[{idx[3:0], 1'b0} + 5'd1]};
    wire signed [15:0] q_i, q_q;
    qpsk_mapper u_qmap (.sym_bits(qbits), .i_sym(q_i), .q_sym(q_q));
    wire [1:0] qhat;
    qpsk_slicer u_qsl (
        .i_sym(q_i + noise_i), .q_sym(q_q + noise_q), .sym_bits(qhat)
    );

    function [1:0] pick_mcs;
        input [2:0] snr_i;
        input [1:0] force_i;
        begin
            if (force_i < 2'd3)
                pick_mcs = force_i;
            else if (snr_i <= 3'd1)
                pick_mcs = 2'd0;
            else if (snr_i <= 3'd3)
                pick_mcs = 2'd1;
            else
                pick_mcs = 2'd2;
        end
    endfunction

    reg        vit_start, vit_valid;
    reg [1:0]  vit_pair;
    wire       vit_busy, vit_done, vit_iv, vit_ib;
    wire [9:0] vit_ix;
    viterbi_hard75 #(.NINFO(32)) u_vit (
        .clk(clk), .rst(rst), .start(vit_start), .rx_valid(vit_valid),
        .rx_pair(vit_pair), .busy(vit_busy), .done(vit_done),
        .info_valid(vit_iv), .info_bit(vit_ib), .info_idx(vit_ix)
    );

    integer ci;

    assign busy     = (state != S_IDLE) && (state != S_DONE);
    assign done     = (state == S_DONE);
    assign mcs_used = mcs_q;
    assign n_ack    = ack_r;
    assign n_fail   = fail_r;
    assign n_uses   = uses_r;
    assign n_good   = good_r;

    always @(posedge clk) begin
        if (rst) begin
            state <= S_IDLE;
            snr_q <= 0;
            mcs_q <= 0;
            max_tx_q <= 1;
            blk <= 0;
            idx <= 0;
            tx_try <= 0;
            saw_low <= 0;
            ack_r <= 0;
            fail_r <= 0;
            uses_r <= 0;
            good_r <= 0;
            crc_rst <= 1;
            crc_en <= 0;
            crc_bit <= 0;
            enc_rst <= 1;
            enc_en <= 0;
            enc_bit <= 0;
            vit_start <= 0;
            vit_valid <= 0;
            vit_pair <= 0;
        end else begin
            crc_rst <= 0;
            crc_en <= 0;
            enc_rst <= 0;
            enc_en <= 0;
            vit_start <= 0;
            vit_valid <= 0;
            if (vit_iv && (vit_ix < 10'd32))
                decb[vit_ix[4:0]] <= vit_ib;

            case (state)
                S_IDLE: begin
                    if (start) begin
                        state    <= S_LOAD;
                        snr_q    <= snr_index;
                        mcs_q    <= pick_mcs(snr_index, mcs_force);
                        max_tx_q <= (max_tx == 0) ? 3'd1 : max_tx;
                        blk      <= 0;
                        ack_r    <= 0;
                        fail_r   <= 0;
                        uses_r   <= 0;
                        good_r   <= 0;
                    end
                end
                S_LOAD: begin
                    state <= S_PAY;
                    idx   <= 0;
                end
                S_PAY: begin
                    pay[idx[3:0]] <= bit_lfsr;
                    if (idx == 7'd15) begin
                        state   <= S_CRC0;
                        crc_rst <= 1;
                    end else
                        idx <= idx + 7'd1;
                end
                // 复位已经在进入本状态的那一拍做完。这里不要再送 pay[0]，
                // 否则第一位会进两次。
                S_CRC0: begin
                    state <= S_CRCIN;
                    idx   <= 0;
                end
                S_CRCIN: begin
                    crc_en  <= 1;
                    crc_bit <= pay[idx[3:0]];
                    if (idx == 7'd15)
                        state <= S_CRCDONE;
                    else
                        idx <= idx + 7'd1;
                end
                S_CRCDONE: begin
                    state <= S_APPEND;
                end
                S_APPEND: begin
                    for (ci = 0; ci < 16; ci = ci + 1) begin
                        word[ci]      <= pay[ci];
                        word[16 + ci] <= crc_reg[15 - ci];
                    end
                    idx <= 0;
                    if (mcs_q == 2'd1) begin
                        state   <= S_ENC0;
                        enc_rst <= 1;
                    end else begin
                        state <= S_TX;
                        tx_try <= 0;
                        for (ci = 0; ci < 68; ci = ci + 1)
                            chase[ci] <= 0;
                    end
                end
                S_ENC0: begin
                    state   <= S_ENC;
                    idx     <= 0;
                    enc_en  <= 1;
                    enc_bit <= word[0];
                end
                S_ENC: begin
                    chan[2 * idx]     <= enc_pair[1];
                    chan[2 * idx + 1] <= enc_pair[0];
                    enc_en <= 1;
                    if (idx == 7'd33) begin
                        state <= S_TX;
                        tx_try <= 0;
                        idx <= 0;
                        for (ci = 0; ci < 68; ci = ci + 1)
                            chase[ci] <= 0;
                    end else begin
                        idx <= idx + 7'd1;
                        enc_bit <= (idx + 7'd1 < 7'd32) ? word[idx + 7'd1] : 1'b0;
                    end
                end
                S_TX: begin
                    uses_r <= uses_r + 32'd1;
                    if (mcs_q == 2'd2) begin
                        chase[{idx[3:0], 1'b0}] <=
                            chase[{idx[3:0], 1'b0}] + (qhat[1] ? -4'sd1 : 4'sd1);
                        chase[{idx[3:0], 1'b0} + 5'd1] <=
                            chase[{idx[3:0], 1'b0} + 5'd1] + (qhat[0] ? -4'sd1 : 4'sd1);
                    end else begin
                        chase[idx] <= chase[idx] + (bpsk_hat ? -4'sd1 : 4'sd1);
                    end
                    if (idx == nsym - 7'd1) begin
                        idx <= 0;
                        if (mcs_q == 2'd1) begin
                            state     <= S_VARM;
                            vit_start <= 1;
                            saw_low   <= 0;
                        end else begin
                            state   <= S_C0;
                            crc_rst <= 1;
                        end
                    end else
                        idx <= idx + 7'd1;
                end
                // start 已在上一拍打出。这一拍 valid 必须保持 0，译码器才不会
                // 把启动沿当成第一对编码比特。34 对在 VFEED 里各送一次。
                S_VARM: begin
                    state <= S_VFEED;
                    idx   <= 0;
                end
                S_VFEED: begin
                    vit_valid <= 1;
                    vit_pair  <= {chase[2 * idx] < 4'sd0, chase[2 * idx + 1] < 4'sd0};
                    if (idx == 7'd33)
                        state <= S_VWAIT;
                    else
                        idx <= idx + 7'd1;
                end
                S_VWAIT: begin
                    if (!vit_done)
                        saw_low <= 1;
                    if (saw_low && vit_done) begin
                        state   <= S_C0;
                        crc_rst <= 1;
                    end
                end
                S_C0: begin
                    state <= S_CCHK;
                    idx   <= 0;
                end
                S_CCHK: begin
                    crc_en <= 1;
                    if (mcs_q == 2'd1)
                        crc_bit <= decb[idx[4:0]];
                    else
                        crc_bit <= (chase[idx] < 4'sd0);
                    if (idx == 7'd31)
                        state <= S_CWAIT;
                    else
                        idx <= idx + 7'd1;
                end
                S_CWAIT: begin
                    state <= S_NEXT;
                end
                S_NEXT: begin
                    if (crc_reg == 16'h0000) begin
                        ack_r  <= ack_r + 16'd1;
                        good_r <= good_r + 32'd16;
                        if (blk == (NBLOCKS - 1))
                            state <= S_DONE;
                        else begin
                            blk   <= blk + 16'd1;
                            idx   <= 0;
                            state <= S_PAY;
                        end
                    end else if (tx_try + 3'd1 < max_tx_q) begin
                        tx_try <= tx_try + 3'd1;
                        idx    <= 0;
                        state  <= S_TX;
                    end else begin
                        fail_r <= fail_r + 16'd1;
                        if (blk == (NBLOCKS - 1))
                            state <= S_DONE;
                        else begin
                            blk   <= blk + 16'd1;
                            idx   <= 0;
                            state <= S_PAY;
                        end
                    end
                end
                S_DONE: begin
                    if (start) begin
                        state    <= S_LOAD;
                        snr_q    <= snr_index;
                        mcs_q    <= pick_mcs(snr_index, mcs_force);
                        max_tx_q <= (max_tx == 0) ? 3'd1 : max_tx;
                        blk      <= 0;
                        ack_r    <= 0;
                        fail_r   <= 0;
                        uses_r   <= 0;
                        good_r   <= 0;
                    end
                end
                default: state <= S_IDLE;
            endcase
        end
    end

endmodule

`default_nettype wire
