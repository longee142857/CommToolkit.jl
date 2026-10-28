`timescale 1ns/1ps
`default_nettype none

// 2×2 复矩阵的 ZF / MMSE。理想 CSI：H 由调用方给出。
//
//   G = H^H H + σ² I
//   x = G^{-1} H^H y
// σ² = noise_var。noise_var = 0 就是 ZF；H 列满秩时它和 MMSE 的 σ²=0 相同。
//
// 刻度：1024 表示 1。QPSK 分量用 724。这组数下 Gram 对角线大约 10^6，
// 2×2 行列式用有符号 64 位。行列式 ≤ 1024 当成奇异，输出 0，singular=1。
// 除法是组合的 64 位有符号除法，商朝 0 截断，再饱和到 16 位。
// 没有流水。H 的实部/虚部绝对值请留在大约 8192 以内，否则 Gram 中间积会溢出。
// 这不是 Julia `pinv` 的逐位结果。
module mimo2x2_detect (
    input  wire signed [15:0] h00_re,
    input  wire signed [15:0] h00_im,
    input  wire signed [15:0] h01_re,
    input  wire signed [15:0] h01_im,
    input  wire signed [15:0] h10_re,
    input  wire signed [15:0] h10_im,
    input  wire signed [15:0] h11_re,
    input  wire signed [15:0] h11_im,
    input  wire signed [31:0] y0_re,
    input  wire signed [31:0] y0_im,
    input  wire signed [31:0] y1_re,
    input  wire signed [31:0] y1_im,
    input  wire [31:0]        noise_var,
    output wire signed [15:0] x0_re,
    output wire signed [15:0] x0_im,
    output wire signed [15:0] x1_re,
    output wire signed [15:0] x1_im,
    output wire               singular
);

    function signed [63:0] mul_hy;
        input signed [15:0] h;
        input signed [31:0] y;
        reg signed [63:0] hh, yy;
        begin
            hh = h;
            yy = y;
            mul_hy = hh * yy;
        end
    endfunction

    function signed [63:0] mul_hh;
        input signed [15:0] a;
        input signed [15:0] b;
        reg signed [63:0] aa, bb;
        begin
            aa = a;
            bb = b;
            mul_hh = aa * bb;
        end
    endfunction

    // z = H^H y。conj(h)*(yr+j yi) 的实部是 hr*yr + hi*yi。
    wire signed [63:0] z0_re =
        mul_hy(h00_re, y0_re) + mul_hy(h00_im, y0_im) +
        mul_hy(h10_re, y1_re) + mul_hy(h10_im, y1_im);
    wire signed [63:0] z0_im =
        mul_hy(h00_re, y0_im) - mul_hy(h00_im, y0_re) +
        mul_hy(h10_re, y1_im) - mul_hy(h10_im, y1_re);
    wire signed [63:0] z1_re =
        mul_hy(h01_re, y0_re) + mul_hy(h01_im, y0_im) +
        mul_hy(h11_re, y1_re) + mul_hy(h11_im, y1_im);
    wire signed [63:0] z1_im =
        mul_hy(h01_re, y0_im) - mul_hy(h01_im, y0_re) +
        mul_hy(h11_re, y1_im) - mul_hy(h11_im, y1_re);

    wire signed [63:0] g00 =
        mul_hh(h00_re, h00_re) + mul_hh(h00_im, h00_im) +
        mul_hh(h10_re, h10_re) + mul_hh(h10_im, h10_im) +
        $signed({32'd0, noise_var});
    wire signed [63:0] g11 =
        mul_hh(h01_re, h01_re) + mul_hh(h01_im, h01_im) +
        mul_hh(h11_re, h11_re) + mul_hh(h11_im, h11_im) +
        $signed({32'd0, noise_var});
    // g01 = conj(h00)*h01 + conj(h10)*h11
    wire signed [63:0] g01_re =
        mul_hh(h00_re, h01_re) + mul_hh(h00_im, h01_im) +
        mul_hh(h10_re, h11_re) + mul_hh(h10_im, h11_im);
    wire signed [63:0] g01_im =
        mul_hh(h00_re, h01_im) - mul_hh(h00_im, h01_re) +
        mul_hh(h10_re, h11_im) - mul_hh(h10_im, h11_re);

    wire signed [63:0] det = g00 * g11 - (g01_re * g01_re + g01_im * g01_im);

    // num = adj(G) z。g11 是实数。
    wire signed [63:0] num0_re = g11 * z0_re - (g01_re * z1_re - g01_im * z1_im);
    wire signed [63:0] num0_im = g11 * z0_im - (g01_re * z1_im + g01_im * z1_re);
    // g10 = conj(g01)
    wire signed [63:0] num1_re = g00 * z1_re - (g01_re * z0_re + g01_im * z0_im);
    wire signed [63:0] num1_im = g00 * z1_im - (g01_re * z0_im - g01_im * z0_re);

    function signed [15:0] sat_div;
        input signed [63:0] num;
        input signed [63:0] den;
        reg signed [63:0] q;
        begin
            if (den == 0)
                q = 0;
            else
                q = num / den;
            if (q > 64'sd32767)
                sat_div = 16'sd32767;
            else if (q < -64'sd32768)
                sat_div = -16'sd32768;
            else
                sat_div = q[15:0];
        end
    endfunction

    assign singular = (det <= 64'sd1024);
    assign x0_re = singular ? 16'sd0 : sat_div(num0_re, det);
    assign x0_im = singular ? 16'sd0 : sat_div(num0_im, det);
    assign x1_re = singular ? 16'sd0 : sat_div(num1_re, det);
    assign x1_im = singular ? 16'sd0 : sat_div(num1_im, det);

endmodule

`default_nettype wire
