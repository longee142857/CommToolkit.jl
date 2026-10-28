// K = 3，八进制 (7, 5)。与 CommToolkit conv_encode 同一套移位：
//   reg = {state[1], state[0], b}
//   g1 = reg[2] xor reg[1] xor reg[0]     // 1 + D + D^2
//   g2 = reg[2] xor reg[0]                // 1 + D^2
//   next = {reg[1], reg[0]}
// 输出对是 {g1, g2}，g1 在前。这个文件要在模块内部 include。

function [1:0] conv75_next;
    input [1:0] s;
    input       b;
    begin
        conv75_next = {s[0], b};
    end
endfunction

function [1:0] conv75_out;
    input [1:0] s;
    input       b;
    reg [2:0] sh;
    begin
        sh = {s, b};
        conv75_out = {^sh, sh[2] ^ sh[0]};
    end
endfunction
