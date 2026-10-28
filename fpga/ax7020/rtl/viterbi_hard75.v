`timescale 1ns/1ps
`default_nettype none

// 终止块上的硬判决 Viterbi。支路度量是一对编码比特的汉明距离。
// 平局用严格小于才替换，状态循环 0→3、每个状态先比特 0 再比特 1，
// 和整数版 viterbi_decode 的书写顺序一致。回溯从路径度量最小的状态开始，
// 平局取编号最小的状态（Julia findmin 的第一个）。不强制结束在状态 0。
// 尾比特两步的判决丢掉，只吐出 NINFO 个信息比特。info_idx 按时间从大到小出现。
// 不是软判决，没有删除符号，所以不做删余。
module viterbi_hard75 #(
    parameter integer NINFO = 64
) (
    input  wire        clk,
    input  wire        rst,
    input  wire        start,
    input  wire        rx_valid,
    input  wire [1:0]  rx_pair,
    output wire        busy,
    output wire        done,
    output wire        info_valid,
    output wire        info_bit,
    output wire [9:0]  info_idx
);

    localparam integer NSYM = NINFO + 2;
    localparam [15:0]  INF  = 16'hFFFE;

    generate
        if ((NINFO < 1) || (NINFO > 256)) begin : g_ninfo
            invalid_viterbi_ninfo_range u_ninfo();
        end
    endgenerate

`include "conv75_trellis.vh"

    localparam [1:0] S_IDLE = 2'd0;
    localparam [1:0] S_RUN  = 2'd1;
    localparam [1:0] S_INIT = 2'd2;
    localparam [1:0] S_TB   = 2'd3;

    reg [1:0]  state;
    reg [9:0]  t_run;
    reg [9:0]  t_tb;
    reg [1:0]  tb_st;
    reg [15:0] pm0, pm1, pm2, pm3;

    // 每步 4 个状态：1 比特判决 + 2 比特前驱。深度 NSYM。
    reg [NSYM*4-1:0] tb_bits;
    reg [NSYM*8-1:0] tb_prev;

    reg [15:0] npm0, npm1, npm2, npm3;
    reg [1:0]  prev0, prev1, prev2, prev3;
    reg        bit0, bit1, bit2, bit3;

    // 不用 task：Icarus 的 always @* 扫不到 task 里读的 rx_pair，支路度量会落后一拍。
    task consider;
        input [1:0]  s;
        input        b;
        input [15:0] pm_s;
        input [1:0]  rx;
        reg [1:0]  ns;
        reg [1:0]  obits;
        reg [15:0] cand;
        begin
            if (pm_s != INF) begin
                ns = conv75_next(s, b);
                obits = conv75_out(s, b);
                cand = pm_s + (rx[1] ^ obits[1]) + (rx[0] ^ obits[0]);
                if ((ns == 2'd0) && (cand < npm0)) begin
                    npm0 = cand; prev0 = s; bit0 = b;
                end
                if ((ns == 2'd1) && (cand < npm1)) begin
                    npm1 = cand; prev1 = s; bit1 = b;
                end
                if ((ns == 2'd2) && (cand < npm2)) begin
                    npm2 = cand; prev2 = s; bit2 = b;
                end
                if ((ns == 2'd3) && (cand < npm3)) begin
                    npm3 = cand; prev3 = s; bit3 = b;
                end
            end
        end
    endtask

    always @* begin
        npm0 = INF; npm1 = INF; npm2 = INF; npm3 = INF;
        prev0 = 2'b00; prev1 = 2'b00; prev2 = 2'b00; prev3 = 2'b00;
        bit0 = 1'b0; bit1 = 1'b0; bit2 = 1'b0; bit3 = 1'b0;
        consider(2'd0, 1'b0, pm0, rx_pair);
        consider(2'd0, 1'b1, pm0, rx_pair);
        consider(2'd1, 1'b0, pm1, rx_pair);
        consider(2'd1, 1'b1, pm1, rx_pair);
        consider(2'd2, 1'b0, pm2, rx_pair);
        consider(2'd2, 1'b1, pm2, rx_pair);
        consider(2'd3, 1'b0, pm3, rx_pair);
        consider(2'd3, 1'b1, pm3, rx_pair);
    end

    function [1:0] best_state;
        input [15:0] a, b, c, d;
        reg [15:0] m;
        begin
            m = a;
            best_state = 2'd0;
            if (b < m) begin m = b; best_state = 2'd1; end
            if (c < m) begin m = c; best_state = 2'd2; end
            if (d < m) begin m = d; best_state = 2'd3; end
        end
    endfunction

    wire [1:0] best = best_state(pm0, pm1, pm2, pm3);
    wire [11:0] slot = (t_tb * 4) + tb_st;

    reg done_r;
    assign busy = (state == S_RUN) || (state == S_INIT) || (state == S_TB);
    assign done = done_r;

    assign info_valid = (state == S_TB) && (t_tb < NINFO);
    assign info_bit   = tb_bits[slot];
    assign info_idx   = t_tb[9:0];

    wire [1:0] prev_now = tb_prev[(slot * 2) +: 2];

    always @(posedge clk) begin
        if (rst) begin
            state  <= S_IDLE;
            done_r <= 1'b0;
            t_run  <= 10'd0;
            t_tb   <= 10'd0;
            tb_st  <= 2'b00;
            pm0 <= INF; pm1 <= INF; pm2 <= INF; pm3 <= INF;
        end else if ((state == S_IDLE || done_r) && start) begin
            state  <= S_RUN;
            done_r <= 1'b0;
            t_run  <= 10'd0;
            pm0 <= 16'd0;
            pm1 <= INF;
            pm2 <= INF;
            pm3 <= INF;
        end else if (state == S_RUN && rx_valid) begin
            pm0 <= npm0; pm1 <= npm1; pm2 <= npm2; pm3 <= npm3;
            tb_bits[(t_run * 4) + 0] <= bit0;
            tb_bits[(t_run * 4) + 1] <= bit1;
            tb_bits[(t_run * 4) + 2] <= bit2;
            tb_bits[(t_run * 4) + 3] <= bit3;
            tb_prev[(((t_run * 4) + 0) * 2) +: 2] <= prev0;
            tb_prev[(((t_run * 4) + 1) * 2) +: 2] <= prev1;
            tb_prev[(((t_run * 4) + 2) * 2) +: 2] <= prev2;
            tb_prev[(((t_run * 4) + 3) * 2) +: 2] <= prev3;
            if (t_run == (NSYM - 1)) begin
                state <= S_INIT;
            end else begin
                t_run <= t_run + 10'd1;
            end
        end else if (state == S_INIT) begin
            tb_st <= best;
            t_tb  <= NSYM - 1;
            state <= S_TB;
        end else if (state == S_TB) begin
            tb_st <= prev_now;
            if (t_tb == 10'd0) begin
                state  <= S_IDLE;
                done_r <= 1'b1;
            end else begin
                t_tb <= t_tb - 9'd1;
            end
        end
    end

endmodule

`default_nettype wire
