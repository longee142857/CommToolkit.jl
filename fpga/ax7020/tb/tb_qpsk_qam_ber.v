`timescale 1ns/1ps
`default_nettype none

// 固定种子下扫 QPSK、16-QAM、64-QAM 的 BER。无噪声必须 0 差错。
// 0 dB 到 8 dB 差错严格下降。区间是趋势检查，不是对闭式公式的验收，
// 也不是和 Julia 浮点仿真比比特。NBITS = 16380，能被 2、4、6 整除。
module tb_qpsk_qam_ber;

    localparam integer NBITS = 16380;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg        rst;
    reg        start_q, start_16, start_64;
    reg [2:0]  snr_q, snr_16, snr_64;
    wire       done_q, done_16, done_64;
    wire [31:0] err_q, err_16, err_64;
    wire [31:0] bit_q, bit_16, bit_64;

    qpsk_qam_awgn_ber #(.MODE(0), .NBITS(NBITS)) u_qpsk (
        .clk(clk), .rst(rst), .start(start_q), .snr_index(snr_q),
        .busy(), .done(done_q), .n_err(err_q), .n_bit(bit_q),
        .tx_bits(), .rx_bits(), .i_tx(), .q_tx(), .i_rx(), .q_rx()
    );
    qpsk_qam_awgn_ber #(.MODE(1), .NBITS(NBITS)) u_qam16 (
        .clk(clk), .rst(rst), .start(start_16), .snr_index(snr_16),
        .busy(), .done(done_16), .n_err(err_16), .n_bit(bit_16),
        .tx_bits(), .rx_bits(), .i_tx(), .q_tx(), .i_rx(), .q_rx()
    );
    qpsk_qam_awgn_ber #(.MODE(2), .NBITS(NBITS)) u_qam64 (
        .clk(clk), .rst(rst), .start(start_64), .snr_index(snr_64),
        .busy(), .done(done_64), .n_err(err_64), .n_bit(bit_64),
        .tx_bits(), .rx_bits(), .i_tx(), .q_tx(), .i_rx(), .q_rx()
    );

    reg [31:0] eq [0:5];
    reg [31:0] e16 [0:5];
    reg [31:0] e64 [0:5];
    integer idx;
    reg [2:0] snr_loop;

    initial begin
        #(NBITS * 40 * 20);
        $display("FAIL: timeout");
        $fatal(1);
    end

    task automatic run_qpsk;
        input [2:0] snr;
        integer db;
        begin
            @(posedge clk);
            snr_q <= snr;
            start_q <= 1'b1;
            @(posedge clk);
            start_q <= 1'b0;
            @(posedge clk);
            if (done_q !== 1'b0) begin
                $display("FAIL: QPSK done stuck after start");
                $fatal(1);
            end
            wait (done_q === 1'b1);
            @(posedge clk);
            if (bit_q !== NBITS) begin
                $display("FAIL: QPSK n_bit %0d", bit_q);
                $fatal(1);
            end
            db = (snr <= 3'd5) ? (snr * 2) : -1;
            if (db < 0)
                $display("QPSK   scale 0          errors %6d / %0d", err_q, bit_q);
            else
                $display("QPSK   Eb/N0 %2d dB     errors %6d / %0d  BER %e", db, err_q, bit_q, err_q * 1.0 / NBITS);
            if (snr <= 3'd5)
                eq[snr] = err_q;
        end
    endtask

    task automatic run_16;
        input [2:0] snr;
        integer db;
        begin
            @(posedge clk);
            snr_16 <= snr;
            start_16 <= 1'b1;
            @(posedge clk);
            start_16 <= 1'b0;
            @(posedge clk);
            if (done_16 !== 1'b0) begin
                $display("FAIL: 16QAM done stuck after start");
                $fatal(1);
            end
            wait (done_16 === 1'b1);
            @(posedge clk);
            if (bit_16 !== NBITS) begin
                $display("FAIL: 16QAM n_bit %0d", bit_16);
                $fatal(1);
            end
            db = (snr <= 3'd5) ? (snr * 2) : -1;
            if (db < 0)
                $display("16QAM  scale 0          errors %6d / %0d", err_16, bit_16);
            else
                $display("16QAM  Eb/N0 %2d dB     errors %6d / %0d  BER %e", db, err_16, bit_16, err_16 * 1.0 / NBITS);
            if (snr <= 3'd5)
                e16[snr] = err_16;
        end
    endtask

    task automatic run_64;
        input [2:0] snr;
        integer db;
        begin
            @(posedge clk);
            snr_64 <= snr;
            start_64 <= 1'b1;
            @(posedge clk);
            start_64 <= 1'b0;
            @(posedge clk);
            if (done_64 !== 1'b0) begin
                $display("FAIL: 64QAM done stuck after start");
                $fatal(1);
            end
            wait (done_64 === 1'b1);
            @(posedge clk);
            if (bit_64 !== NBITS) begin
                $display("FAIL: 64QAM n_bit %0d", bit_64);
                $fatal(1);
            end
            db = (snr <= 3'd5) ? (snr * 2) : -1;
            if (db < 0)
                $display("64QAM  scale 0          errors %6d / %0d", err_64, bit_64);
            else
                $display("64QAM  Eb/N0 %2d dB     errors %6d / %0d  BER %e", db, err_64, bit_64, err_64 * 1.0 / NBITS);
            if (snr <= 3'd5)
                e64[snr] = err_64;
        end
    endtask

    task automatic check_trend;
        input [8*8-1:0] name;
        input [31:0] e0, e1, e2, e3, e4, e5;
        input [31:0] lo0, hi0;
        begin
            if (!(e0 > e1 && e1 > e2 && e2 > e3 && e3 > e4)) begin
                $display("FAIL: %0s not strictly decreasing 0..8 dB (%0d %0d %0d %0d %0d)",
                         name, e0, e1, e2, e3, e4);
                $fatal(1);
            end
            if (e5 > e4) begin
                $display("FAIL: %0s 10 dB above 8 dB (%0d > %0d)", name, e5, e4);
                $fatal(1);
            end
            if (e0 <= lo0 || e0 >= hi0) begin
                $display("FAIL: %0s 0 dB count %0d outside %0d..%0d", name, e0, lo0, hi0);
                $fatal(1);
            end
        end
    endtask

    initial begin
        rst = 1'b1;
        start_q = 1'b0;
        start_16 = 1'b0;
        start_64 = 1'b0;
        snr_q = 3'd0;
        snr_16 = 3'd0;
        snr_64 = 3'd0;
        repeat (4) @(posedge clk);
        rst = 1'b0;
        @(posedge clk);

        run_qpsk(3'd7);
        if (err_q !== 32'd0) begin
            $display("FAIL: noiseless QPSK has %0d errors", err_q);
            $fatal(1);
        end
        run_16(3'd7);
        if (err_16 !== 32'd0) begin
            $display("FAIL: noiseless 16QAM has %0d errors", err_16);
            $fatal(1);
        end
        run_64(3'd7);
        if (err_64 !== 32'd0) begin
            $display("FAIL: noiseless 64QAM has %0d errors", err_64);
            $fatal(1);
        end

        for (idx = 0; idx < 6; idx = idx + 1) begin
            snr_loop = idx;
            run_qpsk(snr_loop);
            run_16(snr_loop);
            run_64(snr_loop);
        end

        // 0 dB：QPSK 收在 4%..15%（教材约 7.9%）；16-QAM 6%..25%；64-QAM 8%..35%。
        check_trend("QPSK", eq[0], eq[1], eq[2], eq[3], eq[4], eq[5],
                    NBITS / 25, (NBITS * 15) / 100);
        check_trend("16QAM", e16[0], e16[1], e16[2], e16[3], e16[4], e16[5],
                    (NBITS * 6) / 100, (NBITS * 25) / 100);
        check_trend("64QAM", e64[0], e64[1], e64[2], e64[3], e64[4], e64[5],
                    (NBITS * 8) / 100, (NBITS * 35) / 100);

        // QPSK 在 4 dB 应明显离开 0 dB 的平台，大约 1%。
        if (eq[2] <= (NBITS * 4) / 1000 || eq[2] >= (NBITS * 4) / 100) begin
            $display("FAIL: QPSK 4 dB count %0d outside 0.4%%..4%%", eq[2]);
            $fatal(1);
        end

        $display("PASS: QPSK, 16-QAM and 64-QAM BER trends");
        $finish;
    end

endmodule

`default_nettype wire
