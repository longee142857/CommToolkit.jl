`timescale 1ns/1ps
`default_nettype none

// 编码器输出和硬 Viterbi 对上几条手算向量，包括 1 个和 2 个编码比特翻转。
// 向量按 conv_encode / 整数 viterbi_decode 的规则算，不从 RTL 反推。
module tb_conv75_unit;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg        enc_rst, enc_en, enc_bit;
    wire [1:0] coded;
    wire [1:0] enc_state;
    conv75_encode u_enc (
        .clk(clk),
        .rst(enc_rst),
        .en(enc_en),
        .bit_in(enc_bit),
        .coded(coded),
        .state(enc_state)
    );

    reg        rst8, start8, valid8;
    reg [1:0]  rx8;
    wire       busy8, done8, iv8, ib8;
    wire [9:0] ix8;
    viterbi_hard75 #(.NINFO(8)) u8 (
        .clk(clk), .rst(rst8), .start(start8), .rx_valid(valid8), .rx_pair(rx8),
        .busy(busy8), .done(done8), .info_valid(iv8), .info_bit(ib8), .info_idx(ix8)
    );

    reg        rst3, start3, valid3;
    reg [1:0]  rx3;
    wire       done3, iv3, ib3;
    wire [9:0] ix3;
    viterbi_hard75 #(.NINFO(3)) u3 (
        .clk(clk), .rst(rst3), .start(start3), .rx_valid(valid3), .rx_pair(rx3),
        .busy(), .done(done3), .info_valid(iv3), .info_bit(ib3), .info_idx(ix3)
    );

    // 信息 1 0 1 1 0 0 1 0，再加两个 0 尾。时间上的第 0 位在最高位。
    localparam [19:0] EXP8 = 20'b1110_0001_0111_1110_1100;
    localparam [7:0]  INFO8 = 8'b1011_0010;
    // 1 0 0 加尾：11 10 11 00 00
    localparam [9:0]  EXP3 = 10'b1110_1100_00;
    localparam [2:0]  INFO3 = 3'b100;

    integer i, ngot;
    reg [7:0] got8;
    reg [2:0] got3;
    reg [19:0] word8;

    task automatic feed8;
        begin
            @(posedge clk);
            start8 <= 1'b1;
            valid8 <= 1'b0;
            @(posedge clk);
            start8 <= 1'b0;
            for (i = 0; i < 10; i = i + 1) begin
                rx8 <= {word8[19 - 2 * i], word8[18 - 2 * i]};
                valid8 <= 1'b1;
                @(posedge clk);
            end
            valid8 <= 1'b0;
            ngot = 0;
            got8 = 8'b0;
            begin : collect8
                integer guard;
                guard = 0;
                while ((done8 !== 1'b1) && (guard < 40)) begin
                    @(posedge clk);
                    #1;
                    if (iv8 === 1'b1) begin
                        got8[ix8] = ib8;
                        ngot = ngot + 1;
                    end
                    guard = guard + 1;
                end
            end
            if (done8 !== 1'b1 || ngot != 8) begin
                $display("FAIL: 8-bit decode done=%b ngot=%0d got=%b", done8, ngot, got8);
                $fatal(1);
            end
        end
    endtask

    initial begin
        enc_rst = 1'b1;
        enc_en = 1'b0;
        enc_bit = 1'b0;
        rst8 = 1'b1;
        start8 = 1'b0;
        valid8 = 1'b0;
        rx8 = 2'b00;
        rst3 = 1'b1;
        start3 = 1'b0;
        valid3 = 1'b0;
        rx3 = 2'b00;
        word8 = EXP8;

        repeat (2) @(posedge clk);
        enc_rst = 1'b0;
        rst8 = 1'b0;
        rst3 = 1'b0;
        @(posedge clk);
        #1;

        // 编码器：状态 0 送 INFO8，再送两个 0。
        for (i = 0; i < 10; i = i + 1) begin
            @(negedge clk);
            if (i < 8)
                enc_bit = INFO8[7 - i];
            else
                enc_bit = 1'b0;
            enc_en = 1'b1;
            #1;
            if (coded !== {EXP8[19 - 2 * i], EXP8[18 - 2 * i]}) begin
                $display("FAIL: encode step %0d got %b state %b bit %b", i, coded, enc_state, enc_bit);
                $fatal(1);
            end
            @(posedge clk);
            @(negedge clk);
            enc_en = 1'b0;
        end
        if (enc_state !== 2'b00) begin
            $display("FAIL: encoder did not return to 0, state %b", enc_state);
            $fatal(1);
        end

        word8 = EXP8;
        feed8;
        // got8[0] 是时间上的第一位。INFO8[7] 才是第一位。
        if (got8 !== 8'b0100_1101) begin
            $display("FAIL: clean decode %b", got8);
            $fatal(1);
        end

        // 0-based 编码比特 2，即 Julia 1-based 下标 3。
        word8 = EXP8;
        word8[17] = ~word8[17];
        feed8;
        if (got8 !== 8'b0100_1101) begin
            $display("FAIL: single flip decode %b", got8);
            $fatal(1);
        end

        // 0-based 9 和 10，即 Julia 1-based 10 和 11。
        word8 = EXP8;
        word8[19 - 9] = ~word8[19 - 9];
        word8[19 - 10] = ~word8[19 - 10];
        feed8;
        if (got8 !== 8'b0100_1101) begin
            $display("FAIL: double flip decode %b", got8);
            $fatal(1);
        end

        // 全 0。
        word8 = 20'b0;
        feed8;
        if (got8 !== 8'b0) begin
            $display("FAIL: zero decode %b", got8);
            $fatal(1);
        end

        // 自由距事件 1 0 0。
        @(posedge clk);
        start3 <= 1'b1;
        @(posedge clk);
        start3 <= 1'b0;
        for (i = 0; i < 5; i = i + 1) begin
            rx3 <= {EXP3[9 - 2 * i], EXP3[8 - 2 * i]};
            valid3 <= 1'b1;
            @(posedge clk);
        end
        valid3 <= 1'b0;
        ngot = 0;
        got3 = 3'b0;
        begin : collect3
            integer guard;
            guard = 0;
            while ((done3 !== 1'b1) && (guard < 20)) begin
                @(posedge clk);
                #1;
                if (iv3 === 1'b1) begin
                    got3[ix3] = ib3;
                    ngot = ngot + 1;
                end
                guard = guard + 1;
            end
        end
        // 时间 0 在 got3[0]，序列 1,0,0 写成 3'b001。
        if (done3 !== 1'b1 || ngot != 3 || got3 !== 3'b001) begin
            $display("FAIL: event 100 decode done=%b n=%0d got=%b", done3, ngot, got3);
            $fatal(1);
        end

        $display("PASS: (7,5) encode and hard Viterbi vectors");
        $finish;
    end

endmodule

`default_nettype wire
