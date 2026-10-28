`default_nettype none

// Fibonacci LFSR。移位约定与 CommToolkit.jl 的 lfsr_sequence 相同：
// 当前输出是 state[0]；反馈是 (state & TAP_MASK) 的奇偶；寄存器朝 bit 0 移动，
// 新比特写入最高位。STEP>1 时在一个时钟内连续推进一步，供噪声源一次取出 8 位。
// 这是教学用线性反馈，不是 ITU-T O.150 测试图案，也不是 NR 扰码。
module lfsr_fib #(
    parameter integer WIDTH    = 15,
    parameter integer TAP_MASK = 15'h4001,
    parameter integer STEP     = 1
) (
    input  wire                  clk,
    input  wire                  rst,
    input  wire                  en,
    input  wire [WIDTH-1:0]      seed,
    output wire                  bit_out,
    output wire [7:0]            sample8
);

    reg [WIDTH-1:0] state;

    // WIDTH >= 8。低 8 位是接下来会移出的 8 个比特，噪声源把它当均匀样本。
    assign bit_out = state[0];
    assign sample8 = state[7:0];

    function [WIDTH-1:0] step_one;
        input [WIDTH-1:0] s;
        reg fb;
        begin
            fb = ^(s & TAP_MASK[WIDTH-1:0]);
            step_one = {fb, s[WIDTH-1:1]};
        end
    endfunction

    reg [WIDTH-1:0] nxt;
    integer si;

    always @* begin
        nxt = state;
        for (si = 0; si < STEP; si = si + 1)
            nxt = step_one(nxt);
    end

    always @(posedge clk) begin
        if (rst)
            state <= seed;
        else if (en)
            state <= nxt;
    end

endmodule

`default_nettype wire
