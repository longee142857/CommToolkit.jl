`timescale 1ns/1ps
`default_nettype none

// 三档 MCS、无噪声 ACK，以及 0 dB 未编码 BPSK 上 max_tx=1 和 max_tx=3 的硬 Chase。
// NBLOCKS=8。计数是教学级，不和 Julia 浮点 simulate_link 比。
module tb_link_mcs;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg        rst, start;
    reg [2:0]  snr;
    reg [1:0]  mcs_force;
    reg [2:0]  max_tx;
    wire       busy, done;
    wire [1:0] mcs_used;
    wire [15:0] n_ack, n_fail;
    wire [31:0] n_uses, n_good;

    link_mcs_harq #(.NBLOCKS(8)) dut (
        .clk(clk),
        .rst(rst),
        .start(start),
        .snr_index(snr),
        .mcs_force(mcs_force),
        .max_tx(max_tx),
        .busy(busy),
        .done(done),
        .mcs_used(mcs_used),
        .n_ack(n_ack),
        .n_fail(n_fail),
        .n_uses(n_uses),
        .n_good(n_good)
    );

    integer guard;
    reg [15:0] fail_once;
    reg [31:0] uses_once;

    task automatic run_link;
        input [2:0] snr_i;
        input [1:0] force_i;
        input [2:0] max_i;
        begin
            @(negedge clk);
            snr       = snr_i;
            mcs_force = force_i;
            max_tx    = max_i;
            start     = 1'b1;
            @(posedge clk);
            @(negedge clk);
            start = 1'b0;
            guard = 0;
            while ((done !== 1'b1) && (guard < 200000)) begin
                @(posedge clk);
                guard = guard + 1;
            end
            if (done !== 1'b1) begin
                $display("FAIL: timeout snr %0d force %0d max %0d", snr_i, force_i, max_i);
                $fatal(1);
            end
            if ((n_ack + n_fail) !== 16'd8) begin
                $display("FAIL: ack+fail %0d+%0d", n_ack, n_fail);
                $fatal(1);
            end
            if (n_good !== (n_ack * 32'd16)) begin
                $display("FAIL: good %0d ack %0d", n_good, n_ack);
                $fatal(1);
            end
        end
    endtask

    initial begin
        rst       = 1'b1;
        start     = 1'b0;
        snr       = 3'd7;
        mcs_force = 2'd0;
        max_tx    = 3'd1;
        repeat (4) @(posedge clk);
        @(negedge clk);
        rst = 1'b0;

        // 无噪声，三档都该整帧通过。符号数：BPSK 32，卷积 68，QPSK 16。
        run_link(3'd7, 2'd0, 3'd1);
        if (mcs_used !== 2'd0 || n_fail !== 0 || n_ack !== 16'd8 || n_uses !== 32'd256) begin
            $display("FAIL: noiseless BPSK mcs %0d ack %0d fail %0d uses %0d",
                     mcs_used, n_ack, n_fail, n_uses);
            $fatal(1);
        end

        run_link(3'd7, 2'd1, 3'd1);
        if (mcs_used !== 2'd1 || n_fail !== 0 || n_ack !== 16'd8 || n_uses !== 32'd544) begin
            $display("FAIL: noiseless conv mcs %0d ack %0d fail %0d uses %0d",
                     mcs_used, n_ack, n_fail, n_uses);
            $fatal(1);
        end

        run_link(3'd7, 2'd2, 3'd1);
        if (mcs_used !== 2'd2 || n_fail !== 0 || n_ack !== 16'd8 || n_uses !== 32'd128) begin
            $display("FAIL: noiseless QPSK mcs %0d ack %0d fail %0d uses %0d",
                     mcs_used, n_ack, n_fail, n_uses);
            $fatal(1);
        end

        // max_tx=0 收成 1，无噪声不应重传。
        run_link(3'd7, 2'd0, 3'd0);
        if (n_fail !== 0 || n_uses !== 32'd256) begin
            $display("FAIL: max_tx 0 uses %0d fail %0d", n_uses, n_fail);
            $fatal(1);
        end

        // 自动档：0/2 dB → 0，4/6 dB → 1，8 dB 和无噪声 → 2。
        run_link(3'd0, 2'd3, 3'd1);
        if (mcs_used !== 2'd0) begin
            $display("FAIL: auto 0 dB picked %0d", mcs_used);
            $fatal(1);
        end
        run_link(3'd1, 2'd3, 3'd1);
        if (mcs_used !== 2'd0) begin
            $display("FAIL: auto 2 dB picked %0d", mcs_used);
            $fatal(1);
        end
        run_link(3'd2, 2'd3, 3'd1);
        if (mcs_used !== 2'd1) begin
            $display("FAIL: auto 4 dB picked %0d", mcs_used);
            $fatal(1);
        end
        run_link(3'd3, 2'd3, 3'd1);
        if (mcs_used !== 2'd1) begin
            $display("FAIL: auto 6 dB picked %0d", mcs_used);
            $fatal(1);
        end
        run_link(3'd4, 2'd3, 3'd1);
        if (mcs_used !== 2'd2) begin
            $display("FAIL: auto 8 dB picked %0d", mcs_used);
            $fatal(1);
        end
        run_link(3'd7, 2'd3, 3'd1);
        if (mcs_used !== 2'd2 || n_fail !== 0 || n_uses !== 32'd128) begin
            $display("FAIL: auto noiseless mcs %0d fail %0d uses %0d",
                     mcs_used, n_fail, n_uses);
            $fatal(1);
        end

        // 0 dB 未编码 BPSK。单次应有失败；三次硬 Chase 的失败数不高于单次，
        // 并且因为重传，信道使用次数更多。
        run_link(3'd0, 2'd0, 3'd1);
        fail_once = n_fail;
        uses_once = n_uses;
        $display("link: Es/N0 0 dB MCS0 max_tx=1 ack=%0d fail=%0d uses=%0d good=%0d",
                 n_ack, n_fail, n_uses, n_good);
        if (n_fail == 0 || n_uses !== 32'd256) begin
            $display("FAIL: 0 dB single-tx fail %0d uses %0d", n_fail, n_uses);
            $fatal(1);
        end

        run_link(3'd0, 2'd0, 3'd3);
        $display("link: Es/N0 0 dB MCS0 max_tx=3 ack=%0d fail=%0d uses=%0d good=%0d",
                 n_ack, n_fail, n_uses, n_good);
        if ((n_fail >= fail_once) || (n_uses <= uses_once)) begin
            $display("FAIL: Chase fail %0d (once %0d) uses %0d (once %0d)",
                     n_fail, fail_once, n_uses, uses_once);
            $fatal(1);
        end

        run_link(3'd0, 2'd1, 3'd1);
        fail_once = n_fail;
        uses_once = n_uses;
        $display("link: Es/N0 0 dB MCS1 max_tx=1 ack=%0d fail=%0d uses=%0d good=%0d",
                 n_ack, n_fail, n_uses, n_good);
        run_link(3'd0, 2'd1, 3'd3);
        $display("link: Es/N0 0 dB MCS1 max_tx=3 ack=%0d fail=%0d uses=%0d good=%0d",
                 n_ack, n_fail, n_uses, n_good);
        if ((n_fail >= fail_once) || (n_uses <= uses_once)) begin
            $display("FAIL: conv Chase fail %0d once %0d uses %0d once %0d",
                     n_fail, fail_once, n_uses, uses_once);
            $fatal(1);
        end

        $display("PASS: MCS/HARQ noiseless, auto select, and hard Chase");
        $finish;
    end

endmodule

`default_nettype wire
