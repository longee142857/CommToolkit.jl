`timescale 1ns/1ps
`default_nettype none

// 直流 1024 的 IFFT 应是 16 个 64，再 FFT 回到 bin 0 的 1024。
// 这锁的是「ifft 除以 N、fft 不除以 N」，不是和 FFTW 的逐位契约。
module tb_fft16;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg rst, start, ifft, wr_en;
    reg [3:0] wr_addr, rd_addr;
    reg signed [15:0] wr_r, wr_i;
    wire busy, done;
    wire signed [15:0] rd_r, rd_i;

    fft16 dut (
        .clk(clk), .rst(rst), .start(start), .ifft(ifft),
        .wr_en(wr_en), .wr_addr(wr_addr), .wr_r(wr_r), .wr_i(wr_i),
        .busy(busy), .done(done), .rd_addr(rd_addr), .rd_r(rd_r), .rd_i(rd_i)
    );

    integer i;
    reg signed [15:0] tr [0:15];
    reg signed [15:0] ti [0:15];

    task automatic run_xform;
        input do_ifft;
        begin
            @(posedge clk);
            ifft <= do_ifft;
            start <= 1'b1;
            @(posedge clk);
            start <= 1'b0;
            ifft <= 1'b0;
            @(posedge clk);
            if (done !== 1'b0) begin
                $display("FAIL: done stayed high after start");
                $fatal(1);
            end
            wait (done === 1'b1);
            @(posedge clk);
        end
    endtask

    initial begin
        rst = 1'b1;
        start = 1'b0;
        ifft = 1'b0;
        wr_en = 1'b0;
        wr_addr = 4'd0;
        wr_r = 16'sd0;
        wr_i = 16'sd0;
        rd_addr = 4'd0;
        repeat (2) @(posedge clk);
        rst = 1'b0;
        @(posedge clk);

        for (i = 0; i < 16; i = i + 1) begin
            @(negedge clk);
            wr_en = 1'b1;
            wr_addr = i[3:0];
            wr_r = (i == 0) ? 16'sd1024 : 16'sd0;
            wr_i = 16'sd0;
            @(posedge clk);
        end
        @(negedge clk);
        wr_en = 1'b0;

        run_xform(1'b1);
        for (i = 0; i < 16; i = i + 1) begin
            rd_addr = i[3:0];
            #1;
            tr[i] = rd_r;
            ti[i] = rd_i;
            if (rd_r !== 16'sd64 || rd_i !== 16'sd0) begin
                $display("FAIL: ifft[%0d] = %0d + j%0d", i, rd_r, rd_i);
                $fatal(1);
            end
        end

        for (i = 0; i < 16; i = i + 1) begin
            @(negedge clk);
            wr_en = 1'b1;
            wr_addr = i[3:0];
            wr_r = tr[i];
            wr_i = ti[i];
            @(posedge clk);
        end
        @(negedge clk);
        wr_en = 1'b0;

        run_xform(1'b0);
        for (i = 0; i < 16; i = i + 1) begin
            rd_addr = i[3:0];
            #1;
            if (i == 0) begin
                if (rd_r !== 16'sd1024 || rd_i !== 16'sd0) begin
                    $display("FAIL: fft dc %0d + j%0d", rd_r, rd_i);
                    $fatal(1);
                end
            end else if (rd_r !== 16'sd0 || rd_i !== 16'sd0) begin
                $display("FAIL: fft bin %0d = %0d + j%0d", i, rd_r, rd_i);
                $fatal(1);
            end
        end

        // 16 个 (+724, +724) 走一圈，硬判决必须仍是 00。
        for (i = 0; i < 16; i = i + 1) begin
            @(negedge clk);
            wr_en = 1'b1;
            wr_addr = i[3:0];
            wr_r = 16'sd724;
            wr_i = 16'sd724;
            @(posedge clk);
        end
        @(negedge clk);
        wr_en = 1'b0;
        run_xform(1'b1);
        for (i = 0; i < 16; i = i + 1) begin
            rd_addr = i[3:0];
            #1;
            tr[i] = rd_r;
            ti[i] = rd_i;
        end
        for (i = 0; i < 16; i = i + 1) begin
            @(negedge clk);
            wr_en = 1'b1;
            wr_addr = i[3:0];
            wr_r = tr[i];
            wr_i = ti[i];
            @(posedge clk);
        end
        @(negedge clk);
        wr_en = 1'b0;
        run_xform(1'b0);
        for (i = 0; i < 16; i = i + 1) begin
            rd_addr = i[3:0];
            #1;
            if (rd_r <= 0 || rd_i <= 0) begin
                $display("FAIL: QPSK bin %0d came back %0d + j%0d", i, rd_r, rd_i);
                $fatal(1);
            end
        end

        $display("PASS: fft16 DC impulse and QPSK roundtrip");
        $finish;
    end

endmodule

`default_nettype wire
