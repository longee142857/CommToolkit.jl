`timescale 1ns/1ps
`default_nettype none

// CRC-16/XMODEM，多项式 0x1021。初值 0，高位先入，不反射，xor-out 0。
// 和 CommToolkit 的 crc16_ccitt 同一套移位：先看 crc[15] xor 输入，再左移，
// 若该比特为 1 则异或 0x1021。ASCII "123456789" 按字节 MSB 先入，余数是 0x31C3。
module crc16_ccitt (
    input  wire        clk,
    input  wire        rst,
    input  wire        en,
    input  wire        bit_in,
    output reg  [15:0] crc
);

    wire mix = crc[15] ^ bit_in;

    always @(posedge clk) begin
        if (rst)
            crc <= 16'h0000;
        else if (en)
            crc <= {crc[14:0], 1'b0} ^ (mix ? 16'h1021 : 16'h0000);
    end

endmodule

`default_nettype wire
