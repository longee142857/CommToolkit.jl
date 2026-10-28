#!/usr/bin/env bash
# 用 Icarus Verilog 跑教学链路。不调用 Vivado，也不需要板子。
# 在 fpga/ax7020 下执行：  bash sim/run_iverilog.sh
# 想看波形：vvp sim/out/tb_bpsk_awgn_ber.vvp +vcd
#   再 gtkwave sim/out/bpsk_awgn_ber.vcd

set -euo pipefail

cd "$(dirname "$0")/.."

if ! command -v iverilog >/dev/null 2>&1 || ! command -v vvp >/dev/null 2>&1; then
    echo "需要 Icarus Verilog（iverilog 和 vvp）。Debian/Ubuntu: sudo apt-get install iverilog" >&2
    exit 1
fi

mkdir -p sim/out

RTL=(
    rtl/lfsr_fib.v
    rtl/bpsk_mapper.v
    rtl/bpsk_slicer.v
    rtl/awgn_clt.v
    rtl/bpsk_awgn_ber.v
    rtl/ax7020_bpsk_top.v
    rtl/qpsk_mapper.v
    rtl/qpsk_slicer.v
    rtl/qam_gray.v
    rtl/qpsk_qam_awgn_ber.v
    rtl/conv75_encode.v
    rtl/viterbi_hard75.v
    rtl/conv_bpsk_awgn_ber.v
    rtl/fft16.v
    rtl/ofdm16_qpsk_ber.v
    rtl/mimo2x2_detect.v
    rtl/mimo2x2_qpsk_ber.v
)

run_one() {
    local tb="$1"
    local name
    name="$(basename "${tb}" .v)"
    local vvp="sim/out/${name}.vvp"
    local log="sim/out/${name}.log"
    echo "== ${name} =="
    iverilog -g2012 -Wall -I rtl -o "${vvp}" "${RTL[@]}" "${tb}"
    set +e
    vvp "${vvp}" > "${log}"
    local code=$?
    set -e
    cat "${log}"
    if [[ ${code} -ne 0 ]]; then
        echo "vvp 退出码 ${code}: ${name}" >&2
        exit 1
    fi
    if grep -q '^FAIL' "${log}"; then
        echo "仿真打印了 FAIL: ${name}" >&2
        exit 1
    fi
    if ! grep -q '^PASS' "${log}"; then
        echo "仿真没有打印 PASS: ${name}" >&2
        exit 1
    fi
}

run_one tb/tb_bpsk_unit.v
run_one tb/tb_bpsk_awgn_ber.v
run_one tb/tb_ax7020_top.v
run_one tb/tb_qpsk_qam_unit.v
run_one tb/tb_qpsk_qam_ber.v
run_one tb/tb_conv75_unit.v
run_one tb/tb_conv_bpsk_ber.v
run_one tb/tb_fft16.v
run_one tb/tb_ofdm16_ber.v
run_one tb/tb_mimo2x2_unit.v
run_one tb/tb_mimo2x2_ber.v

echo "全部 PASS。"
