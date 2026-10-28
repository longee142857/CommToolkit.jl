`timescale 1ns/1ps
`default_nettype none

// 16 点基-2 DIT。输入先做比特反转，输出是自然顺序。
// FFT 的旋转因子是 exp(-j 2π k / m)，不除以 N，对齐 Julia 的 fft。
// IFFT 用共轭旋转因子，并且每一级把蝶形结果算术右移 1 位，四级合起来除以 16，
// 对齐 Julia 的 ifft。因此在没有溢出的时候 fft(ifft(X)) 回到 X，只剩几位舍入。
// 旋转因子用 round(1024 * cos/sin)，乘完再算术右移 10 位。负数朝 -∞。
// 调用方在 start 之前把 16 个样本写进 din；done 之后用 rd_addr 读 mem。
module fft16 (
    input  wire               clk,
    input  wire               rst,
    input  wire               start,
    input  wire               ifft,
    input  wire               wr_en,
    input  wire [3:0]         wr_addr,
    input  wire signed [15:0] wr_r,
    input  wire signed [15:0] wr_i,
    output wire               busy,
    output wire               done,
    input  wire [3:0]         rd_addr,
    output wire signed [15:0] rd_r,
    output wire signed [15:0] rd_i
);

    localparam [1:0] S_IDLE = 2'd0;
    localparam [1:0] S_BREV = 2'd1;
    localparam [1:0] S_RUN  = 2'd2;
    localparam [1:0] S_DONE = 2'd3;

    reg [1:0] state;
    reg       ifft_q;
    reg [3:0] bi;
    reg [4:0] step;

    reg signed [15:0] din_r [0:15];
    reg signed [15:0] din_i [0:15];
    reg signed [15:0] mem_r [0:15];
    reg signed [15:0] mem_i [0:15];

    function [3:0] brev;
        input [3:0] x;
        begin
            brev = {x[0], x[1], x[2], x[3]};
        end
    endfunction

    // FFT 的 W = cos + j sin_fft。IFFT 取共轭，实部相同。
    function signed [15:0] tw_re;
        input [1:0] stage;
        input [2:0] k;
        begin
            case ({stage, k})
                5'b00_000, 5'b01_000, 5'b10_000, 5'b11_000: tw_re = 16'sd1024;
                5'b01_001, 5'b10_010, 5'b11_100: tw_re = 16'sd0;
                5'b10_001, 5'b11_010: tw_re = 16'sd724;
                5'b10_011, 5'b11_110: tw_re = -16'sd724;
                5'b11_001: tw_re = 16'sd946;
                5'b11_011: tw_re = 16'sd392;
                5'b11_101: tw_re = -16'sd392;
                5'b11_111: tw_re = -16'sd946;
                default: tw_re = 16'sd1024;
            endcase
        end
    endfunction

    function signed [15:0] tw_im_fft;
        input [1:0] stage;
        input [2:0] k;
        begin
            case ({stage, k})
                5'b00_000, 5'b01_000, 5'b10_000, 5'b11_000: tw_im_fft = 16'sd0;
                5'b01_001, 5'b10_010, 5'b11_100: tw_im_fft = -16'sd1024;
                5'b10_001, 5'b11_010: tw_im_fft = -16'sd724;
                5'b10_011, 5'b11_110: tw_im_fft = -16'sd724;
                5'b11_001: tw_im_fft = -16'sd392;
                5'b11_011: tw_im_fft = -16'sd946;
                5'b11_101: tw_im_fft = -16'sd946;
                5'b11_111: tw_im_fft = -16'sd392;
                default: tw_im_fft = 16'sd0;
            endcase
        end
    endfunction

    wire [1:0] stage = step[4:3];
    wire [2:0] lin   = step[2:0];
    wire [3:0] mh    = 4'd1 << stage;
    wire [2:0] kidx  = lin & (mh[2:0] - 3'd1);
    wire [2:0] grp   = lin >> stage;
    wire [3:0] base  = (grp * (mh << 1)) + kidx;
    wire [3:0] idx1  = base + mh;

    wire signed [15:0] twr = tw_re(stage, kidx);
    wire signed [15:0] twi = ifft_q ? -tw_im_fft(stage, kidx) : tw_im_fft(stage, kidx);

    wire signed [15:0] ar = mem_r[base];
    wire signed [15:0] ai = mem_i[base];
    wire signed [15:0] br = mem_r[idx1];
    wire signed [15:0] bi_s = mem_i[idx1];

    wire signed [31:0] p_re = (br * twr) - (bi_s * twi);
    wire signed [31:0] p_im = (br * twi) + (bi_s * twr);
    wire signed [15:0] tr = p_re >>> 10;
    wire signed [15:0] ti = p_im >>> 10;

    wire signed [17:0] sum0r = ar + tr;
    wire signed [17:0] sum0i = ai + ti;
    wire signed [17:0] sum1r = ar - tr;
    wire signed [17:0] sum1i = ai - ti;
    wire signed [17:0] y0r_w = ifft_q ? (sum0r >>> 1) : sum0r;
    wire signed [17:0] y0i_w = ifft_q ? (sum0i >>> 1) : sum0i;
    wire signed [17:0] y1r_w = ifft_q ? (sum1r >>> 1) : sum1r;
    wire signed [17:0] y1i_w = ifft_q ? (sum1i >>> 1) : sum1i;

    assign busy = (state == S_BREV) || (state == S_RUN);
    assign done = (state == S_DONE);
    assign rd_r = mem_r[rd_addr];
    assign rd_i = mem_i[rd_addr];

    always @(posedge clk) begin
        if (wr_en && (state == S_IDLE || state == S_DONE)) begin
            din_r[wr_addr] <= wr_r;
            din_i[wr_addr] <= wr_i;
        end
        if (rst) begin
            state  <= S_IDLE;
            ifft_q <= 1'b0;
            bi     <= 4'd0;
            step   <= 5'd0;
        end else if ((state == S_IDLE || state == S_DONE) && start) begin
            state  <= S_BREV;
            ifft_q <= ifft;
            bi     <= 4'd0;
            step   <= 5'd0;
        end else if (state == S_BREV) begin
            mem_r[brev(bi)] <= din_r[bi];
            mem_i[brev(bi)] <= din_i[bi];
            if (bi == 4'd15)
                state <= S_RUN;
            else
                bi <= bi + 4'd1;
        end else if (state == S_RUN) begin
            mem_r[base] <= y0r_w[15:0];
            mem_i[base] <= y0i_w[15:0];
            mem_r[idx1] <= y1r_w[15:0];
            mem_i[idx1] <= y1i_w[15:0];
            if (step == 5'd31)
                state <= S_DONE;
            else
                step <= step + 5'd1;
        end
    end

endmodule

`default_nettype wire
