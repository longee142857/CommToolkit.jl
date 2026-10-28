`timescale 1ns/1ps
`default_nettype none

// 顶层冒烟：复位同步、KEY2 去抖后 SNR 灯加一。不统计 BER（见 tb_bpsk_awgn_ber）。
// 帧长和去抖都缩到仿真能跑完的长度。
module tb_ax7020_top;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg rst_n;
    reg key_n;
    wire [3:0] led;

    ax7020_bpsk_top #(
        .NBITS(64),
        .DEBOUNCE(8),
        .DIV_BIT(3)
    ) dut (
        .clk(clk),
        .rst_n(rst_n),
        .key_n(key_n),
        .led(led)
    );

    initial begin
        rst_n = 1'b1;
        key_n = 1'b1;
        @(posedge clk);
        rst_n = 1'b0;
        repeat (3) @(posedge clk);
        rst_n = 1'b1;

        repeat (120) @(posedge clk);
        #1;
        if (^led === 1'bx) begin
            $display("FAIL: LED is X after reset: %b", led);
            $fatal(1);
        end
        // snr_index = 0 → led[3:1] 灭（高）。led[0] 是心跳。
        if (led[3:1] !== 3'b111) begin
            $display("FAIL: expected SNR 0 on LEDs, got %b", led);
            $fatal(1);
        end

        key_n = 1'b0;
        repeat (30) @(posedge clk);
        #1;
        if (led[1] !== 1'b0 || led[3:2] !== 2'b11) begin
            $display("FAIL: expected SNR 1 after KEY2, got %b", led);
            $fatal(1);
        end

        $display("PASS: AX7020 top reset and SNR key");
        $finish;
    end

endmodule

`default_nettype wire
