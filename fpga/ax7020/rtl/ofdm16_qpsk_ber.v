`timescale 1ns/1ps
`default_nettype none

// 教学 CP-OFDM：16 点 IFFT，循环前缀 4（和 Julia 默认 64/16 同一个 1/4），
// 再 FFT 回来做 QPSK 硬判决。全部 16 个子载波都传数据，包含 DC，没有保护带，
// 也没有导频。这不是 OFDMParams 的默认子载波图。
//
// 平坦信道、符号对齐。不去做 ofdm_sync，也不做 LS/MMSE。H = 1，FFT 后直接判决。
// 时域噪声用刀2 的 QPSK 刻度（k = 2）再算术右移 2 位：未归一化 FFT 把标准差乘 4，
// 先除掉，子载波上的每维噪声才回到 QPSK 的 σ。这不是 apply_channel 的样本功率 SNR。
// snr_index >= 6 时噪声刻度为 0。
module ofdm16_qpsk_ber #(
    parameter integer NSYM = 32
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

    localparam integer NCP = 4;
    localparam integer NSAMP = 16 + NCP;

    localparam [3:0] S_IDLE  = 4'd0;
    localparam [3:0] S_LOAD  = 4'd1;
    localparam [3:0] S_MAP   = 4'd2;
    localparam [3:0] S_I0    = 4'd3;
    localparam [3:0] S_IWAIT = 4'd4;
    localparam [3:0] S_COPY  = 4'd5;
    localparam [3:0] S_TX    = 4'd6;
    localparam [3:0] S_FLOAD = 4'd7;
    localparam [3:0] S_F0    = 4'd8;
    localparam [3:0] S_FWAIT = 4'd9;
    localparam [3:0] S_SLICE = 4'd10;
    localparam [3:0] S_DONE  = 4'd11;

    reg [3:0]  state;
    reg        saw_low;
    reg [2:0]  snr_q;
    reg [15:0] isym;
    reg [4:0]  idx;
    reg [31:0] n_err_r;
    reg [31:0] n_bit_r;
    reg [31:0] tx_bits;

    reg signed [15:0] time_r [0:15];
    reg signed [15:0] time_i [0:15];
    reg signed [15:0] rx_r [0:19];
    reg signed [15:0] rx_i [0:19];

    wire lfsr_rst = rst || (state == S_IDLE) || (state == S_LOAD);
    wire data_en  = (state == S_MAP);
    wire noise_en = (state == S_TX);

    wire data_bit;
    wire [7:0] data_byte;
    lfsr_fib #(
        .WIDTH(15),
        .TAP_MASK(15'h4001),
        .STEP(2)
    ) u_data (
        .clk(clk),
        .rst(lfsr_rst),
        .en(data_en),
        .seed(15'h7FFF),
        .bit_out(data_bit),
        .sample8(data_byte)
    );

    // data_byte[0] 是这一拍最早的比特，作为 QPSK 的 b0。
    wire [1:0] map_bits = {data_byte[0], data_byte[1]};
    wire signed [15:0] map_i, map_q;
    qpsk_mapper u_map (
        .sym_bits(map_bits),
        .i_sym(map_i),
        .q_sym(map_q)
    );

    wire signed [15:0] noise_i, noise_q;
    awgn_clt #(
        .AMP(1024),
        .SEED_SET(0),
        .SCALE_MODE(0)
    ) u_ni (
        .clk(clk),
        .rst(lfsr_rst),
        .en(noise_en),
        .snr_index(snr_q),
        .bits_per_symbol(3'd2),
        .noise(noise_i)
    );
    awgn_clt #(
        .AMP(1024),
        .SEED_SET(1),
        .SCALE_MODE(0)
    ) u_nq (
        .clk(clk),
        .rst(lfsr_rst),
        .en(noise_en),
        .snr_index(snr_q),
        .bits_per_symbol(3'd2),
        .noise(noise_q)
    );

    wire signed [15:0] n_i_s = noise_i >>> 2;
    wire signed [15:0] n_q_s = noise_q >>> 2;

    reg               fft_start;
    reg               fft_ifft;
    reg               fft_wr;
    reg [3:0]         fft_waddr;
    reg signed [15:0] fft_wr_r, fft_wr_i;
    reg [3:0]         fft_raddr;
    wire              fft_busy, fft_done;
    wire signed [15:0] fft_rd_r, fft_rd_i;

    fft16 u_fft (
        .clk(clk),
        .rst(rst),
        .start(fft_start),
        .ifft(fft_ifft),
        .wr_en(fft_wr),
        .wr_addr(fft_waddr),
        .wr_r(fft_wr_r),
        .wr_i(fft_wr_i),
        .busy(fft_busy),
        .done(fft_done),
        .rd_addr(fft_raddr),
        .rd_r(fft_rd_r),
        .rd_i(fft_rd_i)
    );

    wire [1:0] hat_bits;
    qpsk_slicer u_slice (
        .i_sym(fft_rd_r),
        .q_sym(fft_rd_i),
        .sym_bits(hat_bits)
    );

    // CP 取时域末尾 4 点。去 CP 时丢掉前 4 个接收样本。
    wire [4:0] cp_src = (idx < NCP) ? (idx + 5'd12) : (idx - 5'd4);
    wire signed [15:0] tx_r = time_r[cp_src[3:0]] + n_i_s;
    wire signed [15:0] tx_i = time_i[cp_src[3:0]] + n_q_s;

    wire [4:0] body = idx - NCP;
    wire [5:0] bit_at = {idx[3:0], 1'b0};

    assign busy  = (state != S_IDLE) && (state != S_DONE);
    assign done  = (state == S_DONE);
    assign n_err = n_err_r;
    assign n_bit = n_bit_r;

    always @(posedge clk) begin
        if (rst) begin
            state     <= S_IDLE;
            snr_q     <= 3'd0;
            isym      <= 16'd0;
            idx       <= 5'd0;
            n_err_r   <= 32'd0;
            n_bit_r   <= 32'd0;
            fft_start <= 1'b0;
            fft_ifft  <= 1'b0;
            fft_wr    <= 1'b0;
            saw_low   <= 1'b0;
            fft_waddr <= 4'd0;
            fft_wr_r  <= 16'sd0;
            fft_wr_i  <= 16'sd0;
            fft_raddr <= 4'd0;
        end else begin
            fft_start <= 1'b0;
            fft_wr    <= 1'b0;
            case (state)
                S_IDLE: begin
                    if (start) begin
                        state   <= S_LOAD;
                        snr_q   <= snr_index;
                        isym    <= 16'd0;
                        n_err_r <= 32'd0;
                        n_bit_r <= 32'd0;
                    end
                end
                S_LOAD: begin
                    state <= S_MAP;
                    idx   <= 5'd0;
                end
                S_MAP: begin
                    fft_wr    <= 1'b1;
                    fft_waddr <= idx[3:0];
                    fft_wr_r  <= map_i;
                    fft_wr_i  <= map_q;
                    // [2*k] = b0，[2*k+1] = b1。part-select 的高位是拼接的左边。
                    tx_bits[bit_at +: 2] <= {map_bits[0], map_bits[1]};
                    if (idx == 5'd15) begin
                        state <= S_I0;
                    end else begin
                        idx <= idx + 5'd1;
                    end
                end
                S_I0: begin
                    fft_ifft  <= 1'b1;
                    fft_start <= 1'b1;
                    saw_low   <= 1'b0;
                    state     <= S_IWAIT;
                end
                S_IWAIT: begin
                    // done 在 start 被吃掉之前仍可能为 1，先等它掉下来。
                    if (!fft_done)
                        saw_low <= 1'b1;
                    if (saw_low && fft_done) begin
                        state     <= S_COPY;
                        idx       <= 5'd0;
                        fft_raddr <= 4'd0;
                    end
                end
                S_COPY: begin
                    time_r[idx[3:0]] <= fft_rd_r;
                    time_i[idx[3:0]] <= fft_rd_i;
                    if (idx == 5'd15) begin
                        state <= S_TX;
                        idx   <= 5'd0;
                    end else begin
                        idx       <= idx + 5'd1;
                        fft_raddr <= idx[3:0] + 4'd1;
                    end
                end
                S_TX: begin
                    rx_r[idx] <= tx_r;
                    rx_i[idx] <= tx_i;
                    if (idx == (NSAMP - 1)) begin
                        state <= S_FLOAD;
                        idx   <= 5'd0;
                    end else begin
                        idx <= idx + 5'd1;
                    end
                end
                S_FLOAD: begin
                    fft_wr    <= 1'b1;
                    fft_waddr <= idx[3:0];
                    fft_wr_r  <= rx_r[idx + NCP];
                    fft_wr_i  <= rx_i[idx + NCP];
                    if (idx == 5'd15)
                        state <= S_F0;
                    else
                        idx <= idx + 5'd1;
                end
                S_F0: begin
                    fft_ifft  <= 1'b0;
                    fft_start <= 1'b1;
                    saw_low   <= 1'b0;
                    state     <= S_FWAIT;
                end
                S_FWAIT: begin
                    if (!fft_done)
                        saw_low <= 1'b1;
                    if (saw_low && fft_done) begin
                        state     <= S_SLICE;
                        idx       <= 5'd0;
                        fft_raddr <= 4'd0;
                    end
                end
                S_SLICE: begin
                    n_bit_r <= n_bit_r + 32'd2;
                    n_err_r <= n_err_r
                        + {31'd0, hat_bits[1] != tx_bits[{idx[3:0], 1'b0}]}
                        + {31'd0, hat_bits[0] != tx_bits[{idx[3:0], 1'b0} + 5'd1]};
                    if (idx == 5'd15) begin
                        if (isym == (NSYM - 1))
                            state <= S_DONE;
                        else begin
                            isym  <= isym + 16'd1;
                            idx   <= 5'd0;
                            state <= S_MAP;
                        end
                    end else begin
                        idx       <= idx + 5'd1;
                        fft_raddr <= idx[3:0] + 4'd1;
                    end
                end
                S_DONE: begin
                    if (start) begin
                        state   <= S_LOAD;
                        snr_q   <= snr_index;
                        isym    <= 16'd0;
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
