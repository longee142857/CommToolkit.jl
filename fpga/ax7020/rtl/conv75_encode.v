`timescale 1ns/1ps
`default_nettype none

// 终止留给调用方：需要尾比特时把 bit_in 接 0，再 en 两拍。
// 复位后 state = 0。coded 是当前 state 和 bit_in 的组合输出，在 en 的时钟沿之后 state 才前进。
module conv75_encode (
    input  wire       clk,
    input  wire       rst,
    input  wire       en,
    input  wire       bit_in,
    output wire [1:0] coded,
    output reg  [1:0] state
);

`include "conv75_trellis.vh"

    assign coded = conv75_out(state, bit_in);

    always @(posedge clk) begin
        if (rst)
            state <= 2'b00;
        else if (en)
            state <= conv75_next(state, bit_in);
    end

endmodule

`default_nettype wire
