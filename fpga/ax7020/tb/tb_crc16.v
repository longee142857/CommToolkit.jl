`timescale 1ns/1ps
`default_nettype none

// CRC-16/XMODEM。ASCII "123456789" 按字节 MSB 先入，余数锁 0x31C3。
// 再对一条短净荷做追加和整字校验，翻转 1 位后余数必须非 0。
module tb_crc16;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg        rst, en, bit_in;
    wire [15:0] crc;
    crc16_ccitt u (
        .clk(clk),
        .rst(rst),
        .en(en),
        .bit_in(bit_in),
        .crc(crc)
    );

    // "123456789"，最高字节在最高位，每个字节内部也是 MSB 先出。
    localparam [71:0] CHECK = {"1", "2", "3", "4", "5", "6", "7", "8", "9"};
    localparam [7:0]  PAY   = 8'b1011_0010;

    integer i;
    reg [15:0] rem;

    task automatic shift_bit;
        input b;
        begin
            @(negedge clk);
            en     = 1'b1;
            bit_in = b;
            @(posedge clk);
        end
    endtask

    initial begin
        rst    = 1'b1;
        en     = 1'b0;
        bit_in = 1'b0;
        repeat (2) @(posedge clk);
        @(negedge clk);
        rst = 1'b0;

        for (i = 0; i < 72; i = i + 1)
            shift_bit(CHECK[71 - i]);
        @(negedge clk);
        en = 1'b0;
        @(posedge clk);
        #1;
        if (crc !== 16'h31C3) begin
            $display("FAIL: 123456789 remainder %h", crc);
            $fatal(1);
        end

        @(negedge clk);
        rst = 1'b1;
        @(posedge clk);
        @(negedge clk);
        rst = 1'b0;
        for (i = 0; i < 8; i = i + 1)
            shift_bit(PAY[7 - i]);
        @(negedge clk);
        en = 1'b0;
        @(posedge clk);
        #1;
        rem = crc;

        @(negedge clk);
        rst = 1'b1;
        @(posedge clk);
        @(negedge clk);
        rst = 1'b0;
        for (i = 0; i < 8; i = i + 1)
            shift_bit(PAY[7 - i]);
        for (i = 0; i < 16; i = i + 1)
            shift_bit(rem[15 - i]);
        @(negedge clk);
        en = 1'b0;
        @(posedge clk);
        #1;
        if (crc !== 16'h0000) begin
            $display("FAIL: codeword remainder %h", crc);
            $fatal(1);
        end

        @(negedge clk);
        rst = 1'b1;
        @(posedge clk);
        @(negedge clk);
        rst = 1'b0;
        for (i = 0; i < 8; i = i + 1)
            shift_bit((i == 3) ? ~PAY[7 - i] : PAY[7 - i]);
        for (i = 0; i < 16; i = i + 1)
            shift_bit(rem[15 - i]);
        @(negedge clk);
        en = 1'b0;
        @(posedge clk);
        #1;
        if (crc === 16'h0000) begin
            $display("FAIL: one-bit flip was not detected");
            $fatal(1);
        end

        $display("PASS: CRC-16/XMODEM 0x31C3, append, and one-bit flip");
        $finish;
    end

endmodule

`default_nettype wire
