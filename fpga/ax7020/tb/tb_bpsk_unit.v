`timescale 1ns/1ps
`default_nettype none

// 映射、判决、数据 LFSR 前 32 位，以及 AWGN 刻度的定向检查。
// 数据 LFSR 期望值按 lfsr_sequence 的规则：多项式 x^15+x^14+1，种子全 1。
// 噪声期望值是同一套整数公式的前两个样本，第二个是负数，用来锁算术右移。
module tb_bpsk_unit;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg bit_in;
    wire signed [15:0] symbol;
    bpsk_mapper #(.AMP(1024)) u_map (
        .bit_in(bit_in),
        .symbol(symbol)
    );

    reg signed [15:0] slice_in;
    wire              slice_bit;
    bpsk_slicer u_slice (
        .sample(slice_in),
        .bit_out(slice_bit)
    );

    reg        lfsr_rst;
    reg        lfsr_en;
    wire       lfsr_bit;
    lfsr_fib #(
        .WIDTH(15),
        .TAP_MASK(15'h4001),
        .STEP(1)
    ) u_data (
        .clk(clk),
        .rst(lfsr_rst),
        .en(lfsr_en),
        .seed(15'h7FFF),
        .bit_out(lfsr_bit),
        .sample8()
    );

    reg        n_rst;
    reg        n_en;
    wire [7:0] n_byte;
    lfsr_fib #(
        .WIDTH(31),
        .TAP_MASK(31'h1000_0001),
        .STEP(8)
    ) u_noise_lfsr (
        .clk(clk),
        .rst(n_rst),
        .en(n_en),
        .seed(31'h2E4B17C3),
        .bit_out(),
        .sample8(n_byte)
    );

    reg [2:0] snr_index;
    reg [2:0] bits_per_symbol;
    wire signed [15:0] noise;
    awgn_clt #(.AMP(1024), .SEED_SET(0)) u_awgn (
        .clk(clk),
        .rst(n_rst),
        .en(n_en),
        .snr_index(snr_index),
        .bits_per_symbol(bits_per_symbol),
        .noise(noise)
    );

    // 时间上的第 0 位放在 [31]。全 1 种子先移出 15 个 1。
    localparam [31:0] EXP_BITS = 32'b1111_1111_1111_1110_1010_1010_1010_1001;

    integer i;

    task check;
        input integer cond;
        input [8*64-1:0] msg;
        begin
            // !== 1，这样条件为 X 也会失败，不会被当成通过。
            if (cond !== 1) begin
                $display("FAIL: %0s", msg);
                $fatal(1);
            end
        end
    endtask

    initial begin
        bit_in = 1'b0;
        slice_in = 16'sd0;
        lfsr_rst = 1'b1;
        lfsr_en = 1'b0;
        n_rst = 1'b1;
        n_en = 1'b0;
        snr_index = 3'd0;
        bits_per_symbol = 3'd1;

        #1;
        check(symbol === 16'sd1024, "bit 0 maps to +AMP");
        bit_in = 1'b1;
        #1;
        check(symbol === -16'sd1024, "bit 1 maps to -AMP");

        slice_in = 16'sd0;
        #1;
        check(slice_bit === 1'b0, "sample 0 decodes to bit 0");
        slice_in = 16'sd1024;
        #1;
        check(slice_bit === 1'b0, "positive sample decodes to bit 0");
        slice_in = -16'sd1;
        #1;
        check(slice_bit === 1'b1, "sample -1 decodes to bit 1");
        slice_in = -16'sd1024;
        #1;
        check(slice_bit === 1'b1, "sample -AMP decodes to bit 1");

        @(posedge clk);
        #1;
        lfsr_rst = 1'b0;
        lfsr_en = 1'b1;
        n_rst = 1'b0;
        n_en = 1'b1;

        check(lfsr_bit === EXP_BITS[31], "first data bit");
        check(n_byte === 8'hC3, "noise LFSR byte 0");
        check(noise === 16'sd8, "AWGN sample 0 at 0 dB");

        snr_index = 3'd7;
        #1;
        check(noise === 16'sd0, "snr_index 7 is noiseless");
        snr_index = 3'd1;
        #1;
        check(noise === 16'sd6, "AWGN sample 0 at 2 dB");
        snr_index = 3'd0;
        #1;
        check(noise === 16'sd8, "back to 0 dB sample 0");

        @(posedge clk);
        #1;
        check(lfsr_bit === EXP_BITS[30], "second data bit");
        check(n_byte === 8'h17, "noise LFSR byte 1");
        check(noise === -16'sd993, "AWGN sample 1 is negative");

        for (i = 2; i < 32; i = i + 1) begin
            @(posedge clk);
            #1;
            if (lfsr_bit !== EXP_BITS[31-i]) begin
                $display("FAIL: data bit %0d got %b expected %b", i, lfsr_bit, EXP_BITS[31-i]);
                $fatal(1);
            end
        end

        $display("PASS: mapper, slicer, data LFSR, AWGN scale");
        $finish;
    end

endmodule

`default_nettype wire
