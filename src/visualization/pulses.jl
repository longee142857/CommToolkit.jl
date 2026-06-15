# Pulse shaping filters: raised cosine and root-raised cosine

# Raised cosine filter taps
function raised_cosine(beta::Real, span::Int, sps::Int)
    t = range(-span/2, span/2, length=span*sps + 1)
    taps = zeros(length(t))
    for (i, tau) in enumerate(t)
        if tau == 0
            taps[i] = 1.0
        elseif abs(tau) == 1 / (2 * beta)
            taps[i] = (pi / 4) * sinc(1 / (2 * beta))
        else
            val = sinc(tau) * cos(pi * beta * tau) / (1 - (2 * beta * tau)^2)
            taps[i] = val
        end
    end
    taps ./ sum(taps)
end

# Root-raised cosine filter taps
function rrc_filter(beta::Real, span::Int, sps::Int)
    t = range(-span/2, span/2, length=span*sps + 1)
    taps = zeros(length(t))
    for (i, tau) in enumerate(t)
        if tau == 0
            taps[i] = (1 - beta + 4 * beta / pi) / sqrt(sps)
        elseif abs(tau) == 1 / (4 * beta)
            a = beta / sqrt(2)
            b = (1 + 2/pi) * sin(pi / (4 * beta))
            c = (1 - 2/pi) * cos(pi / (4 * beta))
            taps[i] = a * (b + c) / sqrt(sps)
        else
            num = sin(pi * tau * (1 - beta)) + 4 * beta * tau * cos(pi * tau * (1 + beta))
            den = pi * tau * (1 - (4 * beta * tau)^2)
            taps[i] = den != 0 ? num / den / sqrt(sps) : 0.0
        end
    end
    taps ./ sqrt(sum(abs2, taps))
end

# Upsample and apply pulse shaping filter
function pulse_shape(symbols::AbstractVector{<:Number}, taps::AbstractVector{<:Real}, sps::Int)
    up = zeros(eltype(symbols), length(symbols) * sps)
    up[1:sps:end] .= symbols
    # Manual convolution
    n, m = length(up), length(taps)
    result = zeros(eltype(up), n + m - 1)
    for i in 1:n, j in 1:m
        result[i + j - 1] += up[i] * taps[j]
    end
    result[1:n]
end

# Sinc function
sinc(x::Real) = x == 0 ? 1.0 : sin(pi * x) / (pi * x)
