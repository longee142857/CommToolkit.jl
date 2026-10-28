`default_nettype none

// AX7020 PL 顶层封装。只把教学链路接到板载 50 MHz、两个按键和四只 LED。
// 没有 Block Design，没有 PS、AXI、DDR 或 UART，也没有 ILA。
// 误码计数留在 n_err / n_bit，不上引脚。本文件没有在板上跑过。
//
// KEY1（低有效）复位。KEY2（低有效）在 0,2,4,6,8,10 dB 之间循环。
// 测量自动一帧接一帧。LED 低电平点亮：
//   led[0] 心跳（时钟在走）
//   led[1] snr_index[0]
//   led[2] snr_index[1]
//   led[3] snr_index[2]
// 索引 0..5 就是上面六个 Eb/N0。两只按键都要先同步。
module ax7020_bpsk_top #(
    parameter integer NBITS     = 8_388_608,
    parameter integer DEBOUNCE  = 1_000_000,
    parameter integer DIV_BIT   = 25
) (
    input  wire       clk,
    input  wire       rst_n,
    input  wire       key_n,
    output wire [3:0] led
);

    reg [1:0] rst_sync;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            rst_sync <= 2'b00;
        else
            rst_sync <= {rst_sync[0], 1'b1};
    end
    wire rst = ~rst_sync[1];

    reg        key_sync0, key_sync1, key_stable, key_stable_d;
    reg [31:0] db_cnt;

    always @(posedge clk) begin
        if (rst) begin
            key_sync0   <= 1'b1;
            key_sync1   <= 1'b1;
            key_stable  <= 1'b1;
            key_stable_d <= 1'b1;
            db_cnt      <= 32'd0;
        end else begin
            key_sync0    <= key_n;
            key_sync1    <= key_sync0;
            key_stable_d <= key_stable;
            if (key_sync1 == key_stable) begin
                db_cnt <= 32'd0;
            end else if (db_cnt == (DEBOUNCE - 1)) begin
                key_stable <= key_sync1;
                db_cnt     <= 32'd0;
            end else begin
                db_cnt <= db_cnt + 32'd1;
            end
        end
    end

    wire key_press = key_stable_d & ~key_stable;

    reg [2:0]  snr_index;
    reg        start;
    reg        arm;
    reg [25:0] hb;
    wire       done;

    always @(posedge clk) begin
        if (rst) begin
            snr_index <= 3'd0;
            start     <= 1'b0;
            arm       <= 1'b1;
            hb        <= 26'd0;
        end else begin
            hb <= hb + 26'd1;
            if (key_press) begin
                if (snr_index == 3'd5)
                    snr_index <= 3'd0;
                else
                    snr_index <= snr_index + 3'd1;
            end
            if (arm) begin
                start <= 1'b1;
                arm   <= 1'b0;
            end else begin
                start <= done;
            end
        end
    end

    bpsk_awgn_ber #(
        .NBITS(NBITS),
        .AMP(1024)
    ) u_link (
        .clk(clk),
        .rst(rst),
        .start(start),
        .snr_index(snr_index),
        .busy(),
        .done(done),
        .n_err(),
        .n_bit(),
        .tx_bit(),
        .rx_bit(),
        .symbol(),
        .noise(),
        .sample()
    );

    wire heartbeat = hb[DIV_BIT];
    assign led = {~snr_index, ~heartbeat};

endmodule

`default_nettype wire
