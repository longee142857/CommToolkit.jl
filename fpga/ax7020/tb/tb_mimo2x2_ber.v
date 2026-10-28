`timescale 1ns/1ps
`default_nettype none

// 固定耦合 2×2 上的 QPSK。ZF 和 MMSE 都要求无噪声为 0，并且 0→4 dB 下降。
// 不要求 MMSE 每一档都低于 ZF：1024 比特太短。
module tb_mimo2x2_ber;

    localparam integer NVEC = 256;
    localparam integer NTOTAL = NVEC * 4;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg rst, start_z, start_m;
    reg [2:0] snr_z, snr_m;
    wire done_z, done_m;
    wire [31:0] err_z, err_m, bit_z, bit_m;

    mimo2x2_qpsk_ber #(.MODE(0), .NVEC(NVEC)) u_zf (
        .clk(clk), .rst(rst), .start(start_z), .snr_index(snr_z),
        .busy(), .done(done_z), .n_err(err_z), .n_bit(bit_z)
    );
    mimo2x2_qpsk_ber #(.MODE(1), .NVEC(NVEC)) u_mmse (
        .clk(clk), .rst(rst), .start(start_m), .snr_index(snr_m),
        .busy(), .done(done_m), .n_err(err_m), .n_bit(bit_m)
    );

    reg [31:0] ez [0:4];
    reg [31:0] em [0:4];
    integer k;
    reg [2:0] snr_loop;

    initial begin
        #(NTOTAL * 30 * 20);
        $display("FAIL: timeout");
        $fatal(1);
    end

    task automatic run_zf;
        input [2:0] snr;
        integer db;
        begin
            @(posedge clk);
            snr_z <= snr;
            start_z <= 1'b1;
            @(posedge clk);
            start_z <= 1'b0;
            @(posedge clk);
            if (done_z !== 1'b0) begin
                $display("FAIL: ZF done stuck");
                $fatal(1);
            end
            wait (done_z === 1'b1);
            @(posedge clk);
            if (bit_z !== NTOTAL) begin
                $display("FAIL: ZF n_bit %0d", bit_z);
                $fatal(1);
            end
            db = (snr <= 3'd4) ? snr * 2 : -1;
            if (db < 0)
                $display("ZF   scale 0       errors %6d / %0d", err_z, bit_z);
            else
                $display("ZF   Eb/N0 %2d dB  errors %6d / %0d  BER %e", db, err_z, bit_z, err_z * 1.0 / NTOTAL);
            if (snr <= 3'd4)
                ez[snr] = err_z;
        end
    endtask

    task automatic run_mmse;
        input [2:0] snr;
        integer db;
        begin
            @(posedge clk);
            snr_m <= snr;
            start_m <= 1'b1;
            @(posedge clk);
            start_m <= 1'b0;
            @(posedge clk);
            if (done_m !== 1'b0) begin
                $display("FAIL: MMSE done stuck");
                $fatal(1);
            end
            wait (done_m === 1'b1);
            @(posedge clk);
            if (bit_m !== NTOTAL) begin
                $display("FAIL: MMSE n_bit %0d", bit_m);
                $fatal(1);
            end
            db = (snr <= 3'd4) ? snr * 2 : -1;
            if (db < 0)
                $display("MMSE scale 0       errors %6d / %0d", err_m, bit_m);
            else
                $display("MMSE Eb/N0 %2d dB  errors %6d / %0d  BER %e", db, err_m, bit_m, err_m * 1.0 / NTOTAL);
            if (snr <= 3'd4)
                em[snr] = err_m;
        end
    endtask

    task automatic check_trend;
        input [8*8-1:0] name;
        input [31:0] e0, e1, e2, e3, e4;
        begin
            if (!(e0 > e1 && e1 > e2)) begin
                $display("FAIL: %0s not decreasing 0>2>4 (%0d %0d %0d)", name, e0, e1, e2);
                $fatal(1);
            end
            if (e3 > e2 || e4 > e3) begin
                $display("FAIL: %0s high SNR rose (%0d %0d %0d)", name, e2, e3, e4);
                $fatal(1);
            end
            if (e0 <= (NTOTAL / 50) || e0 >= ((NTOTAL * 25) / 100)) begin
                $display("FAIL: %0s 0 dB %0d outside 2%%..25%%", name, e0);
                $fatal(1);
            end
        end
    endtask

    initial begin
        rst = 1'b1;
        start_z = 1'b0;
        start_m = 1'b0;
        snr_z = 3'd0;
        snr_m = 3'd0;
        repeat (4) @(posedge clk);
        rst = 1'b0;
        @(posedge clk);

        run_zf(3'd7);
        if (err_z !== 0) begin
            $display("FAIL: noiseless ZF %0d", err_z);
            $fatal(1);
        end
        run_mmse(3'd7);
        if (err_m !== 0) begin
            $display("FAIL: noiseless MMSE %0d", err_m);
            $fatal(1);
        end

        for (k = 0; k < 5; k = k + 1) begin
            snr_loop = k;
            run_zf(snr_loop);
            run_mmse(snr_loop);
        end

        check_trend("ZF", ez[0], ez[1], ez[2], ez[3], ez[4]);
        check_trend("MMSE", em[0], em[1], em[2], em[3], em[4]);

        $display("PASS: 2x2 ZF and MMSE BER trends");
        $finish;
    end

endmodule

`default_nettype wire
