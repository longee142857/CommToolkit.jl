# CommToolkit.jl

A Julia communication systems simulation toolkit for learning and research.

## Features

- **Modulation**: BPSK, QPSK, 16QAM with Gray coding
- **Channel models**: AWGN, flat Rayleigh fading, TDL multipath
- **Channel coding**: rate-1/2 convolutional code with Viterbi decoder
- **Equalization**: ZF and MMSE linear equalizers for frequency-selective channels
- **OFDM**: complete transceiver chain with LS/MMSE channel estimation and STO/CFO synchronization

## Project Structure

```
CommToolkit.jl/
├── src/
│   ├── CommToolkit.jl       # Module entry point
│   ├── modems/              # BPSK, QPSK, 16QAM
│   ├── channels/            # AWGN, Rayleigh, TDL
│   ├── coding/              # Convolutional code + Viterbi
│   ├── equalizers/          # ZF, MMSE
│   ├── ofdm/                # OFDM tx/rx/sync
│   └── utils/               # BER, SNR conversion
├── demos/                   # Runnable simulation scripts
│   ├── ber_bpsk.jl          # BPSK over AWGN
│   ├── ber_qpsk.jl          # QPSK over AWGN
│   ├── ber_16qam.jl         # 16QAM over AWGN
│   ├── ber_bpsk_rayleigh.jl # BPSK over Rayleigh fading
│   ├── ber_conv_coding.jl   # Convolutional coding gain
│   ├── ber_equalizer.jl     # ZF vs MMSE comparison
│   ├── ofdm_ber.jl          # OFDM with channel estimation
│   └── ber_comparison.jl    # All schemes compared
├── test/runtests.jl         # 22 unit tests
└── Project.toml
```

## Quick Start

```julia
using Pkg; Pkg.add("https://github.com/longee142857/CommToolkit.jl")
using CommToolkit

# Simulate BPSK over AWGN
bits = rand(0:1, 10000)
syms = modulate(bits, BPSK())
rx = awgn(syms, 10, 1)  # 10 dB Eb/N0
rx_bits = demodulate(rx, BPSK())
error_rate = ber(bits, rx_bits)
```

## Dependencies

- Julia 1.10+
- FFTW (OFDM only)
- Plots.jl (demos only)

## License

MIT
