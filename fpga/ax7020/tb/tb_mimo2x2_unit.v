`timescale 1ns/1ps
`default_nettype none

// 手算过的 2x2 向量。正交矩阵的 ZF 回到 ±724。
// 耦合矩阵的 ZF 同样回到符号，MMSE（noise_var = 524288）幅度缩小但符号不变。
module tb_mimo2x2_unit;

    reg signed [15:0] h00r, h00i, h01r, h01i, h10r, h10i, h11r, h11i;
    reg signed [31:0] y0r, y0i, y1r, y1i;
    reg [31:0] nvar;
    wire signed [15:0] x0r, x0i, x1r, x1i;
    wire singular;

    mimo2x2_detect dut (
        .h00_re(h00r), .h00_im(h00i), .h01_re(h01r), .h01_im(h01i),
        .h10_re(h10r), .h10_im(h10i), .h11_re(h11r), .h11_im(h11i),
        .y0_re(y0r), .y0_im(y0i), .y1_re(y1r), .y1_im(y1i),
        .noise_var(nvar),
        .x0_re(x0r), .x0_im(x0i), .x1_re(x1r), .x1_im(x1i),
        .singular(singular)
    );

    initial begin
        h00i = 0;
        h01i = 0;
        h10i = 0;
        h11i = 0;
        y0i = 0;
        y1i = 0;
        nvar = 0;

        // 正交：列 (724,724) 与 (724,-724)。x=(724,724) 时 y=(1048352, 0)。
        h00r = 16'sd724;
        h01r = 16'sd724;
        h10r = 16'sd724;
        h11r = -16'sd724;
        y0r = 32'sd1048352;
        y1r = 0;
        #1;
        if (singular !== 1'b0 || x0r !== 16'sd724 || x0i !== 0 || x1r !== 16'sd724 || x1i !== 0) begin
            $display("FAIL: orth ++ (%0d,%0d) (%0d,%0d) sing %b", x0r, x0i, x1r, x1i, singular);
            $fatal(1);
        end

        // x=(-724,-724),(724,-724)。y=(0,-1048352),(-1048352,0)。
        y0r = 0;
        y0i = -32'sd1048352;
        y1r = -32'sd1048352;
        y1i = 0;
        #1;
        if (singular !== 1'b0 || x0r !== -16'sd724 || x0i !== -16'sd724 ||
            x1r !== 16'sd724 || x1i !== -16'sd724) begin
            $display("FAIL: orth qpsk (%0d,%0d) (%0d,%0d)", x0r, x0i, x1r, x1i);
            $fatal(1);
        end

        // 耦合实矩阵。y 对应 x=(724, -724)。
        h00r = 16'sd1024;
        h01r = 16'sd256;
        h10r = 16'sd128;
        h11r = 16'sd1024;
        y0r = 32'sd556032;
        y0i = 0;
        y1r = -32'sd648704;
        y1i = 0;
        #1;
        if (singular !== 1'b0 || x0r !== 16'sd724 || x0i !== 0 || x1r !== -16'sd724 || x1i !== 0) begin
            $display("FAIL: coupled ZF (%0d,%0d) (%0d,%0d) sing %b", x0r, x0i, x1r, x1i, singular);
            $fatal(1);
        end

        nvar = 32'd524288;
        #1;
        if (singular !== 1'b0 || x0r !== 16'sd409 || x0i !== 0 || x1r !== -16'sd416 || x1i !== 0) begin
            $display("FAIL: coupled MMSE (%0d,%0d) (%0d,%0d) sing %b", x0r, x0i, x1r, x1i, singular);
            $fatal(1);
        end

        // 两列相同，det = 0。
        h00r = 16'sd1024;
        h01r = 16'sd1024;
        h10r = 16'sd1024;
        h11r = 16'sd1024;
        y0r = 32'sd1;
        y1r = 32'sd1;
        nvar = 0;
        #1;
        if (singular !== 1'b1 || x0r !== 0 || x1r !== 0) begin
            $display("FAIL: singular (%0d,%0d) sing %b", x0r, x1r, singular);
            $fatal(1);
        end

        $display("PASS: 2x2 ZF/MMSE fixed vectors");
        $finish;
    end

endmodule

`default_nettype wire
