`timescale 1ns/1ps
`default_nettype none

// 刀2 定向检查：QPSK / 16-QAM / 64-QAM 的几个绝对坐标、全部星座点的
// 硬判决往返，以及 awgn_clt 在 k = 2/4/6 时的第一拍刻度。
// 坐标按 round(1024/sqrt(...)) 手算，不从 RTL 里抄函数。
module tb_qpsk_qam_unit;

    reg clk;
    initial clk = 1'b0;
    always #5 clk = ~clk;

    reg [1:0] qpsk_bits;
    wire signed [15:0] qpsk_i, qpsk_q;
    wire [1:0] qpsk_hat;
    qpsk_mapper u_qpsk_map (
        .sym_bits(qpsk_bits),
        .i_sym(qpsk_i),
        .q_sym(qpsk_q)
    );
    qpsk_slicer u_qpsk_slice (
        .i_sym(qpsk_i),
        .q_sym(qpsk_q),
        .sym_bits(qpsk_hat)
    );

    reg [5:0] qam_bits;
    wire signed [15:0] i16, q16, i64, q64;
    wire [5:0] hat16, hat64;
    qam_mapper #(.M(16)) u_map16 (
        .sym_bits(qam_bits),
        .i_sym(i16),
        .q_sym(q16)
    );
    qam_slicer #(.M(16)) u_slice16 (
        .i_sym(i16),
        .q_sym(q16),
        .sym_bits(hat16)
    );
    qam_mapper #(.M(64)) u_map64 (
        .sym_bits(qam_bits),
        .i_sym(i64),
        .q_sym(q64)
    );
    qam_slicer #(.M(64)) u_slice64 (
        .i_sym(i64),
        .q_sym(q64),
        .sym_bits(hat64)
    );

    reg        n_rst;
    reg        n_en;
    reg [2:0]  snr_index;
    reg [2:0]  bps;
    wire signed [15:0] n0, n1;
    awgn_clt #(.AMP(1024), .SEED_SET(0)) u_n0 (
        .clk(clk),
        .rst(n_rst),
        .en(n_en),
        .snr_index(snr_index),
        .bits_per_symbol(bps),
        .noise(n0)
    );
    awgn_clt #(.AMP(1024), .SEED_SET(1)) u_n1 (
        .clk(clk),
        .rst(n_rst),
        .en(n_en),
        .snr_index(snr_index),
        .bits_per_symbol(bps),
        .noise(n1)
    );

    integer i;

    task expect_qpsk;
        input [1:0] bits;
        input signed [15:0] ei;
        input signed [15:0] eq;
        begin
            qpsk_bits = bits;
            #1;
            if (qpsk_i !== ei || qpsk_q !== eq || qpsk_hat !== bits) begin
                $display("FAIL: QPSK %b -> (%0d,%0d) hat %b", bits, qpsk_i, qpsk_q, qpsk_hat);
                $fatal(1);
            end
        end
    endtask

    task expect_qam16;
        input [3:0] bits;
        input signed [15:0] ei;
        input signed [15:0] eq;
        begin
            qam_bits = {2'b00, bits};
            #1;
            if (i16 !== ei || q16 !== eq || hat16[3:0] !== bits) begin
                $display("FAIL: 16QAM %b -> (%0d,%0d) hat %b", bits, i16, q16, hat16);
                $fatal(1);
            end
        end
    endtask

    task expect_qam64;
        input [5:0] bits;
        input signed [15:0] ei;
        input signed [15:0] eq;
        begin
            qam_bits = bits;
            #1;
            if (i64 !== ei || q64 !== eq || hat64 !== bits) begin
                $display("FAIL: 64QAM %b -> (%0d,%0d) hat %b", bits, i64, q64, hat64);
                $fatal(1);
            end
        end
    endtask

    initial begin
        qpsk_bits = 2'b00;
        qam_bits = 6'b000000;
        n_rst = 1'b1;
        n_en = 1'b0;
        snr_index = 3'd0;
        bps = 3'd2;

        // b0 在高位。00 在第一象限，11 在第三象限。A = 724。
        expect_qpsk(2'b00, 16'sd724, 16'sd724);
        expect_qpsk(2'b01, 16'sd724, -16'sd724);
        expect_qpsk(2'b10, -16'sd724, 16'sd724);
        expect_qpsk(2'b11, -16'sd724, -16'sd724);

        // 16-QAM，U = 324。高两位是 I 的 Gray。
        expect_qam16(4'b0000, -16'sd972, -16'sd972);
        expect_qam16(4'b0100, -16'sd324, -16'sd972);
        expect_qam16(4'b0010, -16'sd972, 16'sd972);
        expect_qam16(4'b1001, 16'sd972, -16'sd324);
        expect_qam16(4'b1111, 16'sd324, 16'sd324);

        expect_qam64(6'b000000, -16'sd1106, -16'sd1106);
        expect_qam64(6'b100000, 16'sd1106, -16'sd1106);
        expect_qam64(6'b000100, -16'sd1106, 16'sd1106);
        expect_qam64(6'b111111, 16'sd474, 16'sd474);
        expect_qam64(6'b010101, -16'sd158, 16'sd790);

        for (i = 0; i < 4; i = i + 1) begin
            qpsk_bits = i[1:0];
            #1;
            if (qpsk_hat !== qpsk_bits) begin
                $display("FAIL: QPSK roundtrip %b -> %b", qpsk_bits, qpsk_hat);
                $fatal(1);
            end
        end
        for (i = 0; i < 16; i = i + 1) begin
            qam_bits = i[5:0];
            #1;
            if (hat16[3:0] !== qam_bits[3:0]) begin
                $display("FAIL: 16QAM roundtrip %b -> %b", qam_bits[3:0], hat16);
                $fatal(1);
            end
        end
        for (i = 0; i < 64; i = i + 1) begin
            qam_bits = i[5:0];
            #1;
            if (hat64 !== qam_bits) begin
                $display("FAIL: 64QAM roundtrip %b -> %b", qam_bits, hat64);
                $fatal(1);
            end
        end

        @(posedge clk);
        #1;
        n_rst = 1'b0;
        // 种子组 0 的第一拍中心化结果是 +3。k = 2 时 MULT = 512，噪声 = 6。
        // 种子组 1 的中心化结果是 -660，噪声 = -1320。两组不能相同。
        if (n0 !== 16'sd6 || n1 !== -16'sd1320) begin
            $display("FAIL: QPSK noise samples %0d %0d", n0, n1);
            $fatal(1);
        end
        bps = 3'd4;
        #1;
        if (n0 !== 16'sd4) begin
            $display("FAIL: 16QAM scale sample %0d", n0);
            $fatal(1);
        end
        bps = 3'd6;
        #1;
        if (n0 !== 16'sd3) begin
            $display("FAIL: 64QAM scale sample %0d", n0);
            $fatal(1);
        end
        bps = 3'd1;
        #1;
        if (n0 !== 16'sd8) begin
            $display("FAIL: BPSK scale regressed to %0d", n0);
            $fatal(1);
        end

        $display("PASS: QPSK, 16-QAM, 64-QAM map/slice and noise scales");
        $finish;
    end

endmodule

`default_nettype wire
