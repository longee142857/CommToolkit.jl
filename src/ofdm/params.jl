# OFDM system parameters

struct OFDMParams
    N_fft::Int           # FFT size
    N_data::Int          # number of data subcarriers
    N_pilot::Int         # number of pilot subcarriers
    N_cp::Int            # cyclic prefix length
    mod_order::Int       # modulation order (4=QPSK, 16=16QAM)
    bits_per_symbol::Int # bits per symbol
    pilot_spacing::Int   # pilot spacing in subcarriers
    pilot_indices::Vector{Int}
    data_indices::Vector{Int}
    null_indices::Vector{Int}
end

function OFDMParams(; N_fft::Int=64, N_cp::Int=16, mod_order::Int=4,
                      n_guard::Int=6, pilot_spacing::Int=6)
    active_indices = collect(n_guard+1 : N_fft - n_guard)
    dc_idx = N_fft ÷ 2 + 1
    used_indices = filter(x -> x != dc_idx, active_indices)

    pilot_positions = collect(pilot_spacing:pilot_spacing:length(used_indices))
    pilot_indices = sort(used_indices[pilot_positions])

    data_indices = sort(setdiff(used_indices, pilot_indices))
    bits_per_sym = Int(log2(mod_order))

    OFDMParams(N_fft, length(data_indices), length(pilot_indices),
               N_cp, mod_order, bits_per_sym, pilot_spacing,
               pilot_indices, data_indices, [dc_idx])
end

# Derive modulation scheme from parameters
modulation(params::OFDMParams) = params.mod_order == 4 ? QPSK() : QAM(params.mod_order)

# Simple LCG for reproducible random numbers (no Random dependency)
function _lcg_rand(seed::Int)
    r = (seed * 1103515245 + 12345) & 0x7fffffff
    Float64(r) / 2147483648.0
end
