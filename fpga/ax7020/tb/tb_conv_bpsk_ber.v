`timescale 1ns/1ps
`default_nettype none

// (7,5) + BPSK + AWGN + 硬 Viterbi 的信息比特 BER。无噪声必须为 0。
// 0 dB 到 4 dB 差错下降。不是和 Julia 浮点或联合界的数值契约。
module tb_conv_bpsk_ber;

    localparam integer NINFO   = 64;
    localparam integer NBLOCKS = 32;
    localparam integer NTOTAL  = NINFO * NBLOCKS;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg        rst, start;
    reg [2:0]  snr_index;
    wire       busy, done;
    wire [31:0] n_err, n_bit;

    conv_bpsk_awgn_ber #(
        .NINFO(NINFO),
        .NBLOCKS(NBLOCKS)
    ) dut (
        .clk(clk),
        .rst(rst),
        .start(start),
        .snr_index(snr_index),
        .busy(busy),
        .done(done),
        .n_err(n_err),
        .n_bit(n_bit)
    );

    reg [31:0] errs [0:4];
    integer idx;
    reg [2:0] snr_loop;

    initial begin
        #(NTOTAL * 40 * 20);
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
                $display("FAIL: done stuck after start");
                $fatal(1);
            end
            wait (done === 1'b1);
            @(posedge clk);
            if (n_bit !== NTOTAL) begin
                $display("FAIL: n_bit %0d expected %0d", n_bit, NTOTAL);
                $fatal(1);
            end
            db = (snr <= 3'd4) ? (snr * 2) : -1;
            if (db < 0)
                $display("conv75 scale 0       errors %6d / %0d", n_err, n_bit);
            else
                $display("conv75 Eb/N0 %2d dB  errors %6d / %0d  BER %e",
                         db, n_err, n_bit, n_err * 1.0 / NTOTAL);
            if (snr <= 3'd4)
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
            $display("FAIL: noiseless coded link has %0d info errors", n_err);
            $fatal(1);
        end

        for (idx = 0; idx < 5; idx = idx + 1) begin
            snr_loop = idx;
            run_point(snr_loop);
        end

        if (!(errs[0] > errs[1] && errs[1] > errs[2])) begin
            $display("FAIL: coded BER not decreasing 0>2>4 dB (%0d %0d %0d)",
                     errs[0], errs[1], errs[2]);
            $fatal(1);
        end
        if (errs[3] > errs[2] || errs[4] > errs[3]) begin
            $display("FAIL: high SNR rose (%0d %0d %0d)", errs[2], errs[3], errs[4]);
            $fatal(1);
        end
        // 0 dB 信息比特差错应离开 0，也远不到一半。2048 比特的 2%..40%。
        if (errs[0] <= (NTOTAL / 50) || errs[0] >= ((NTOTAL * 40) / 100)) begin
            $display("FAIL: 0 dB count %0d outside 2%%..40%%", errs[0]);
            $fatal(1);
        end

        $display("PASS: coded BPSK hard-Viterbi BER trend");
        $finish;
    end

endmodule

`default_nettype wire
