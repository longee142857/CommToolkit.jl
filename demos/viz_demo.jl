import CommToolkit: modulate, BPSK, QPSK, QAM, awgn
import CommToolkit: raised_cosine, rrc_filter, pulse_shape
import CommToolkit: plot_constellation, plot_eye_diagram, plot_pulse_response
import CommToolkit: snr_db_to_linear

# === 1. Pulse shapes ===
rc = raised_cosine(0.35, 8, 8)
rrc = rrc_filter(0.35, 8, 8)
p_pulse = plot_pulse_response(rc, title="Raised Cosine (β=0.35)")
Plots.savefig(joinpath(@__DIR__, "viz_pulse_rc.png"))

p_rrc = plot_pulse_response(rrc, title="RRC (β=0.35)")
Plots.savefig(joinpath(@__DIR__, "viz_pulse_rrc.png"))

# === 2. Eye diagram ===
sps = 8
syms = modulate(rand(0:1, 500), BPSK())
shaped = pulse_shape(syms, rrc, sps)
p_eye = plot_eye_diagram(shaped, sps, title="BPSK Eye Diagram (RRC, β=0.35)")
Plots.savefig(joinpath(@__DIR__, "viz_eye.png"))

# === 3. Constellation diagrams ===
# BPSK
tx = rand(0:1, 1000)
syms = modulate(tx, BPSK())
rx = awgn(syms, 12, 1)
p1 = plot_constellation(rx, BPSK(), title="BPSK at 12 dB")
Plots.savefig(joinpath(@__DIR__, "viz_constellation_bpsk.png"))

# QPSK
syms = modulate(rand(0:1, 2000), QPSK())
rx = awgn(syms, 12, 2)
p2 = plot_constellation(rx, QPSK(), title="QPSK at 12 dB")
Plots.savefig(joinpath(@__DIR__, "viz_constellation_qpsk.png"))

# 16QAM
syms = modulate(rand(0:1, 4000), QAM(16))
rx = awgn(syms, 16, 4)
p3 = plot_constellation(rx, QAM(16), title="16QAM at 16 dB")
Plots.savefig(joinpath(@__DIR__, "viz_constellation_16qam.png"))

println("Visualization demo complete.")
