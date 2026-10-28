`default_nettype none

// BPSK 硬判决，概念上对应 demodulate：负数判为 1，0 和正数判为 0。
// 只看实部（本链路的样本本来就是实数）。没有软信息，没有载波恢复。
module bpsk_slicer (
    input  wire signed [15:0] sample,
    output wire               bit_out
);

    assign bit_out = (sample < 16'sd0);

endmodule

`default_nettype wire
