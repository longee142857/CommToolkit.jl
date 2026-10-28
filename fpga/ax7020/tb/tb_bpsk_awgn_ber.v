`timescale 1ns/1ps
`default_nettype none

// 固定种子下扫 Eb/N0。无噪声必须 0 差错；0 dB 到 8 dB 差错严格下降；
// 10 dB 不超过 8 dB。区间是宽松的趋势检查，不是对 Q 函数的验收。
module tb_bpsk_awgn_ber;

    localparam integer NBITS = 16384;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg        rst;
    reg        start;
    reg [2:0]  snr_index;
    wire       busy;
    wire       done;
    wire [31:0] n_err;
    wire [31:0] n_bit;
    wire       tx_bit;
    wire       rx_bit;
    wire signed [15:0] symbol;
    wire signed [15:0] noise;
    wire signed [15:0] sample;

    bpsk_awgn_ber #(
        .NBITS(NBITS),
        .AMP(1024)
    ) dut (
        .clk(clk),
        .rst(rst),
        .start(start),
        .snr_index(snr_index),
        .busy(busy),
        .done(done),
        .n_err(n_err),
        .n_bit(n_bit),
        .tx_bit(tx_bit),
        .rx_bit(rx_bit),
        .symbol(symbol),
        .noise(noise),
        .sample(sample)
    );

    reg [31:0] errs [0:5];
    reg [2:0]  snr_loop;
    integer idx;

    initial begin
        if ($test$plusargs("vcd")) begin
            $dumpfile("sim/out/bpsk_awgn_ber.vcd");
            $dumpvars(0, tb_bpsk_awgn_ber);
        end
    end

    initial begin
        #(NBITS * 20 * 20);
        $display("FAIL: timeout");
        $fatal(1);
    end

    task automatic run_point;
        input [2:0] snr;
        integer db;
        begin
            @(posedge clk);
            snr_index <= snr;
            start <= 1'b1;
            @(posedge clk);
            start <= 1'b0;
            @(posedge clk);
            if (done !== 1'b0) begin
                $display("FAIL: done did not clear after start (snr_index %0d)", snr);
                $fatal(1);
            end
            wait (done === 1'b1);
            @(posedge clk);
            if (n_bit !== NBITS) begin
                $display("FAIL: n_bit %0d, expected %0d", n_bit, NBITS);
                $fatal(1);
            end
            db = (snr <= 3'd5) ? (snr * 2) : -1;
            if (db < 0)
                $display("scale 0 (noiseless)  errors %6d / %0d", n_err, n_bit);
            else
                $display("Eb/N0 %2d dB          errors %6d / %0d  BER %e", db, n_err, n_bit, n_err * 1.0 / NBITS);
            if (snr <= 3'd5)
                errs[snr] = n_err;
        end
    endtask

    initial begin
        rst = 1'b1;
        start = 1'b0;
        snr_index = 3'd0;
        repeat (4) @(posedge clk);
        rst = 1'b0;
        @(posedge clk);

        run_point(3'd7);
        if (n_err !== 32'd0) begin
            $display("FAIL: noiseless run has %0d errors", n_err);
            $fatal(1);
        end

        for (idx = 0; idx < 6; idx = idx + 1) begin
            snr_loop = idx;
            run_point(snr_loop);
        end

        if (!(errs[0] > errs[1] && errs[1] > errs[2] && errs[2] > errs[3] && errs[3] > errs[4])) begin
            $display("FAIL: error counts are not strictly decreasing from 0 dB to 8 dB");
            $fatal(1);
        end
        if (errs[5] > errs[4]) begin
            $display("FAIL: 10 dB has more errors than 8 dB (%0d > %0d)", errs[5], errs[4]);
            $fatal(1);
        end
        // 16384 比特：0 dB 理论约 7.9%，这里收在 5%..12%；4 dB 理论约 1.25%，收在 0.4%..4%。
        if (errs[0] <= (NBITS / 20) || errs[0] >= ((NBITS * 12) / 100)) begin
            $display("FAIL: 0 dB error count %0d is outside 5%%..12%%", errs[0]);
            $fatal(1);
        end
        if (errs[2] <= ((NBITS * 4) / 1000) || errs[2] >= ((NBITS * 4) / 100)) begin
            $display("FAIL: 4 dB error count %0d is outside 0.4%%..4%%", errs[2]);
            $fatal(1);
        end

        $display("PASS: noiseless is 0; BER falls from 0 dB through 8 dB");
        $finish;
    end

endmodule

`default_nettype wire
