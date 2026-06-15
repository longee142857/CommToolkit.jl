# Convolutional code: encoder + Viterbi decoder
#
# Rate 1/2, constraint length K=3 (4 states)
# Generators: G1 = 7 (111₂), G2 = 5 (101₂)

const CONV_G1 = 0b111  # polynomial: 1 + D + D^2
const CONV_G2 = 0b101  # polynomial: 1 + D^2
const CONV_K = 3        # constraint length
const CONV_N_STATES = 4 # 2^(K-1) states

# Precompute next state and output for each (state, input_bit)
function _conv_trellis()
    next_state = zeros(Int, CONV_N_STATES, 2)
    output = zeros(Int, CONV_N_STATES, 2)

    for s in 0:CONV_N_STATES-1
        for b in 0:1
            reg = (s << 1) | b  # shift in new bit
            # Generator outputs
            g1 = count_ones(reg & CONV_G1) & 1
            g2 = count_ones(reg & CONV_G2) & 1
            next_state[s+1, b+1] = reg & (CONV_N_STATES - 1)
            output[s+1, b+1] = (g1 << 1) | g2
        end
    end
    return next_state, output
end

const CONV_NEXT, CONV_OUT = _conv_trellis()

# Encode bits with convolutional code
function conv_encode(bits::AbstractVector{<:Integer})
    n = length(bits)
    coded = Vector{Int}(undef, 2 * (n + CONV_K - 1))  # + tail bits
    state = 0

    # Encode data bits
    idx = 1
    for b in bits
        out = CONV_OUT[state+1, b+1]
        coded[idx]   = (out >> 1) & 1
        coded[idx+1] = out & 1
        state = CONV_NEXT[state+1, b+1]
        idx += 2
    end

    # Tail bits (flush to state 0)
    for _ in 1:CONV_K-1
        out = CONV_OUT[state+1, 1]  # input 0
        coded[idx]   = (out >> 1) & 1
        coded[idx+1] = out & 1
        state = CONV_NEXT[state+1, 1]
        idx += 2
    end

    return coded
end

# Viterbi decoder: hard-decision
function viterbi_decode(coded::AbstractVector{<:Integer})
    n_coded = length(coded)
    @assert n_coded % 2 == 0 "coded length must be even"
    n_sym = n_coded ÷ 2

    # Trellis: path metric, traceback
    n_states = CONV_N_STATES
    INF = typemax(Int) ÷ 2
    pm = fill(INF, n_states)
    pm[1] = 0  # start at state 0
    tb = zeros(Int, n_sym, n_states)

    # Precompute branch metrics for each state/input
    for t in 1:n_sym
        rx0 = coded[2t-1]
        rx1 = coded[2t]
        new_pm = fill(INF, n_states)

        for s in 0:n_states-1
            if pm[s+1] == INF
                continue
            end
            for b in 0:1
                out = CONV_OUT[s+1, b+1]
                expected0 = (out >> 1) & 1
                expected1 = out & 1
                # Hamming distance
                bm = (rx0 ⊻ expected0) + (rx1 ⊻ expected1)
                ns = CONV_NEXT[s+1, b+1]
                cand = pm[s+1] + bm
                if cand < new_pm[ns+1]
                    new_pm[ns+1] = cand
                    tb[t, ns+1] = s
                end
            end
        end
        pm = new_pm
    end

    # Traceback (from best final state)
    _, best_state = findmin(pm)
    best_state -= 1

    # Remove tail bits from length
    n_bits = n_sym - (CONV_K - 1)
    bits = Vector{Int}(undef, n_bits)

    for t in n_sym:-1:1
        if t <= n_bits
            # Determine input bit from prev_state → state transition
            prev = tb[t, best_state+1]
            # Which input from prev would reach best_state?
            # Check both inputs to prev_state
            for b in 0:1
                if CONV_NEXT[prev+1, b+1] == best_state
                    bits[t] = b
                    break
                end
            end
        end
        best_state = tb[t, best_state+1]
    end

    return bits
end
