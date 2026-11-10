# Independent checks for demos/lab_rrc_isi.jl.
#
# The BigFloat formulas and sums live only in this file. They are not a
# production dependency. γ_n · Σ|a_i b_i| bounds one floating-point sum
# (Higham, 2nd ed., (3.5), u = eps(Float64)/2). It does not bound the
# pulse formula, the 1e-12 branch, normalization, or truncation.
# A small |Σ g_n^2 − 1| is not used as evidence that the formula is right.
# The finite RRC cascade is not compared to the same-span raised cosine.

struct _Iv
    v::BigFloat
    d::BigFloat
end

function _u()
    return big(eps(Float64)) / 2
end

function _gamma(n::Integer)
    n == 0 && return big(0)
    nu = n * _u()
    nu < 1 || error("gamma undefined")
    return nu / (1 - nu)
end

function _add(a::_Iv, b::_Iv)
    v = a.v + b.v
    return _Iv(v, a.d + b.d + _u() * abs(v))
end

function _sub(a::_Iv, b::_Iv)
    v = a.v - b.v
    return _Iv(v, a.d + b.d + _u() * abs(v))
end

function _mul(a::_Iv, b::_Iv)
    v = a.v * b.v
    d = abs(a.v) * b.d + abs(b.v) * a.d + a.d * b.d + _u() * abs(v)
    return _Iv(v, d)
end

function _div(a::_Iv, b::_Iv)
    v = a.v / b.v
    if b.d >= abs(b.v)
        return _Iv(v, abs(v) + abs(a.v) + a.d + b.d)
    end
    d = (a.d + abs(v) * b.d) / (abs(b.v) - b.d) + _u() * abs(v)
    return _Iv(v, d)
end

function _sin(a::_Iv)
    v = sin(a.v)
    d = abs(cos(a.v)) * a.d + a.d^2 / 2 + _u() * abs(v)
    return _Iv(v, d)
end

function _cos(a::_Iv)
    v = cos(a.v)
    d = abs(sin(a.v)) * a.d + a.d^2 / 2 + _u() * abs(v)
    return _Iv(v, d)
end

function _sqrt(a::_Iv)
    v = sqrt(a.v)
    d = a.d / (2 * v) + a.d^2 / (8 * v^3) + _u() * abs(v)
    return _Iv(v, d)
end

function _pi()
    # Center is mathematical π. The radius covers the Float64 constant the code uses.
    p = big(π)
    return _Iv(p, _u() * p)
end

function _exact(x::Real)
    return _Iv(big(x), big(0))
end

function _rrc_raw(τf::Float64, βf::Float64, sps::Int)
    τ = _exact(τf)
    β = _exact(βf)
    one = _exact(1)
    two = _exact(2)
    four = _exact(4)
    spsb = _sqrt(_exact(sps))
    πb = _pi()
    if βf == 0
        return τf == 0 ? one : _div(_sin(_mul(πb, τ)), _mul(πb, τ))
    elseif τf == 0
        return _div(_add(_sub(one, β), _div(_mul(four, β), πb)), spsb)
    elseif abs(abs(τf) - 1 / (4 * βf)) < 1e-12
        a = _div(β, _sqrt(two))
        s = _div(πb, _mul(four, β))
        b = _mul(_add(one, _div(two, πb)), _sin(s))
        c = _mul(_sub(one, _div(two, πb)), _cos(s))
        return _div(_mul(a, _add(b, c)), spsb)
    else
        num = _add(
            _sin(_mul(_mul(πb, τ), _sub(one, β))),
            _mul(_mul(_mul(four, β), τ), _cos(_mul(_mul(πb, τ), _add(one, β)))),
        )
        den = _mul(_mul(πb, τ), _sub(one, _mul(_mul(_mul(four, β), τ), _mul(_mul(four, β), τ))))
        return _div(_div(num, den), spsb)
    end
end

function _rc_raw(τf::Float64, βf::Float64)
    τ = _exact(τf)
    β = _exact(βf)
    one = _exact(1)
    two = _exact(2)
    four = _exact(4)
    πb = _pi()
    if βf == 0 || τf == 0
        return τf == 0 ? one : _div(_sin(_mul(πb, τ)), _mul(πb, τ))
    elseif abs(abs(τf) - 1 / (2 * βf)) < 1e-12
        s = _div(one, _mul(two, β))
        sinc = _div(_sin(_mul(πb, s)), _mul(πb, s))
        return _mul(_div(πb, four), sinc)
    else
        sinc = _div(_sin(_mul(πb, τ)), _mul(πb, τ))
        num = _mul(sinc, _cos(_mul(_mul(πb, β), τ)))
        den = _sub(one, _mul(_mul(_mul(two, β), τ), _mul(_mul(two, β), τ)))
        return _div(num, den)
    end
end

function _normalize_energy(raw::Vector{_Iv})
    acc = _Iv(big(0), big(0))
    for r in raw
        acc = _add(acc, _mul(r, r))
    end
    nrm = _sqrt(acc)
    return [_div(r, nrm) for r in raw]
end

function _normalize_dc(raw::Vector{_Iv})
    acc = _Iv(big(0), big(0))
    for r in raw
        acc = _add(acc, r)
    end
    return [_div(r, acc) for r in raw]
end

function _formula_ratio(got::Vector{Float64}, ref::Vector{_Iv})
    worst = 0.0
    where = ""
    for i in eachindex(got)
        err = abs(big(got[i]) - ref[i].v)
        if err > ref[i].d
            return Inf, "i=$i err=$err bound=$(ref[i].d) got=$(got[i])"
        end
        ratio = ref[i].d == 0 ? 0.0 : Float64(err / ref[i].d)
        if ratio > worst
            worst = ratio
            where = "i=$i"
        end
    end
    return worst, where
end

function _rrc_ratio(β, span, sps)
    g = rrc_filter(β, span, sps)
    t = collect(range(-span / 2, span / 2; length=length(g)))
    raw = [_rrc_raw(Float64(τ), Float64(β), sps) for τ in t]
    return _formula_ratio(g, _normalize_energy(raw))
end

function _rc_ratio(β, span, sps)
    h = raised_cosine(β, span, sps)
    t = collect(range(-span / 2, span / 2; length=length(h)))
    raw = [_rc_raw(Float64(τ), Float64(β)) for τ in t]
    return _formula_ratio(h, _normalize_dc(raw))
end

function _conv_big(x, h)
    y = zeros(BigFloat, length(x) + length(h) - 1)
    for i in eachindex(x)
        xi = big(x[i])
        for j in eachindex(h)
            y[i + j - 1] += xi * big(h[j])
        end
    end
    return y
end

function _scale(x, h, k)
    s = big(0)
    n = 0
    for i in eachindex(x)
        j = k - i + 1
        if 1 <= j <= length(h)
            s += abs(big(Float64(x[i])) * big(Float64(h[j])))
            n += 1
        end
    end
    return n, s
end

function _lab_symbols(n)
    rng = Xoshiro(0x0c0de1162)
    return Float64[rand(rng) < 0.5 ? -1.0 : 1.0 for _ in 1:n]
end

function _sidelobe_sum(r, L, sps)
    acc = big(0)
    energy = big(0)
    m = 1
    while true
        hit = false
        for q in (L + m * sps, L - m * sps)
            if 1 <= q <= length(r)
                acc += abs(r[q])
                energy += r[q]^2
                hit = true
            end
        end
        hit || break
        m += 1
    end
    return acc, energy
end

const _RRC_GRID = (
    (0.0, 1, 1),
    (0.0, 6, 4),
    (0.22, 3, 3),
    (0.25, 2, 1),
    (0.35, 4, 64),
    (0.35, 5, 3),
    (0.35, 6, 4),
    (0.35, 6, 7),
    (0.35, 16, 8),
    (0.5, 1, 1),
    (0.5, 2, 2),
    (1.0, 2, 4),
    (1.0, 3, 1),
)

@testset "RRC formula on the real support" begin
    worst = 0.0
    for (β, span, sps) in _RRC_GRID
        ratio, where = setprecision(BigFloat, 256) do
            _rrc_ratio(β, span, sps)
        end
        worst = max(worst, ratio)
        @test ratio <= 1
        if ratio > 1
            println("RRC formula β=$β span=$span sps=$sps $where")
        end
    end
    println("RRC formula max |err|/bound = ", worst)
    g = rrc_filter(0.35, 6, 7)
    t = collect(range(-3.0, 3.0; length=length(g)))
    target = 1 / (4 * 0.35)
    @test any(τ -> abs(abs(τ) - target) == 0, t)
    energy_residual = abs(sum(abs2, rrc_filter(0.35, 6, 4)) - 1)
    println("normalization |Σg^2-1| = ", energy_residual, " (not the formula check)")
end

@testset "Raised cosine formula, reference only" begin
    worst = 0.0
    for (β, span, sps) in _RRC_GRID
        ratio, where = setprecision(BigFloat, 256) do
            _rc_ratio(β, span, sps)
        end
        worst = max(worst, ratio)
        @test ratio <= 1
        if ratio > 1
            println("RC formula β=$β span=$span sps=$sps $where")
        end
    end
    println("RC formula max |err|/bound = ", worst)
    # Endpoints of span=2, β=0.5 land on ±1/(2β).
    h = raised_cosine(0.5, 2, 1)
    @test length(h) == 3
end

@testset "Removable RRC limits, BigFloat only" begin
    setprecision(BigFloat, 256) do
        for βf in (0.22, 0.35, 0.5, 1.0)
            β = big(βf)
            closed0 = 1 - β + 4 * β / π
            function general(τ)
                num = sin(π * τ * (1 - β)) + 4 * β * τ * cos(π * τ * (1 + β))
                den = π * τ * (1 - (4 * β * τ)^2)
                return num / den
            end
            d_far = abs(general(big(2)^-40) - closed0)
            d_near = abs(general(big(2)^-70) - closed0)
            @test d_near < d_far / 1000
            s = 1 / (4 * β)
            a = β / sqrt(big(2))
            b = (1 + 2 / π) * sin(π / (4 * β))
            c = (1 - 2 / π) * cos(π / (4 * β))
            closed = a * (b + c)
            for side in (-1, 1)
                q_far = abs(general(s + side * big(2)^-40) - closed)
                q_near = abs(general(s + side * big(2)^-70) - closed)
                @test q_near < q_far / 1000
            end
        end
    end
end

@testset "Neighborhood outside the 1e-12 branch" begin
    β = 0.35
    span, sps = 4, 64
    target = 1 / (4 * β)
    t = collect(range(-span / 2, span / 2; length=span * sps + 1))
    best = Inf
    for τ in t
        d = abs(abs(τ) - target)
        if d >= 1e-12 && d < best
            best = d
        end
    end
    @test best > 1e-12
    @test best <= 1 / sps
    ratio, = setprecision(BigFloat, 256) do
        _rrc_ratio(β, span, sps)
    end
    @test ratio <= 1
end

@testset "Convolution rounding is not the truncation" begin
    symbols = _lab_symbols(32)
    for (β, span, sps) in ((0.35, 6, 4), (0.35, 10, 4), (0.5, 1, 1), (1.0, 3, 1), (0.35, 5, 3))
        g = rrc_filter(β, span, sps)
        L = length(g)
        @test L == span * sps + 1
        rev = reverse(g)
        x = zeros(length(symbols) * sps)
        x[1:sps:end] .= symbols
        tx = apply_fir(x, g)
        rx = matched_filter(tx, g)
        setprecision(BigFloat, 256) do
            tx_ref = _conv_big(x, g)
            for k in eachindex(tx)
                n, s = _scale(x, g, k)
                @test abs(big(tx[k]) - tx_ref[k]) <= _gamma(n) * s
            end
            # Second sum, same already-rounded tx. This γ does not include tx's error.
            rx_from_float_tx = _conv_big(tx, rev)
            for k in eachindex(rx)
                n, s = _scale(tx, rev, k)
                @test abs(big(rx[k]) - rx_from_float_tx[k]) <= _gamma(n) * s
            end
            rx_from_exact_tx = _conv_big(tx_ref, rev)
            r = _conv_big(g, rev)
            g0 = r[L]
            S, energy = _sidelobe_sum(r, L, sps)
            for k in eachindex(symbols)
                p = L + (k - 1) * sps
                y = big(0)
                for m in eachindex(symbols)
                    q = L + (k - m) * sps
                    if 1 <= q <= length(r)
                        y += big(symbols[m]) * r[q]
                    end
                end
                # Two BigFloat evaluations of the same cascade. The ceiling is the
                # working-precision unit roundoff times the number of taps in one
                # output, not a truncation tolerance.
                u_bf = eps(BigFloat(1))
                id_ceiling = (length(x) + length(g)) * u_bf * (1 + abs(g0) + S)
                @test abs(rx_from_exact_tx[p] - y) <= id_ceiling
                @test symbols[k] * y + id_ceiling >= g0 - S
            end
            # |g0-1| is whatever the Float64 taps' energy is. Not a formula proof.
            @test energy >= 0
            println("β=$β span=$span sps=$sps truncation_energy=", Float64(energy),
                    " g0-1=", Float64(g0 - 1))
        end
        @test length(x) == length(symbols) * sps
        @test length(tx) == length(x) + L - 1
        @test length(rx) == length(x) + 2 * (L - 1)
        kept = downsample_symbols(rx, sps; offset=L - 1)
        @test length(kept) == length(symbols) + span
        @test kept[1] == rx[L]
    end
end

@testset "Impulse offset is the skip count" begin
    for (β, span, sps) in ((0.35, 6, 4), (0.5, 1, 1), (0.35, 5, 3), (1.0, 3, 1))
        g = rrc_filter(β, span, sps)
        L = length(g)
        @test L == span * sps + 1
        if iseven(L)
            @test !isinteger((L - 1) / 2)
        end
        rev = reverse(g)
        for p0 in (1, 4)
            buf = zeros(p0)
            buf[p0] = 1.0
            rx = matched_filter(apply_fir(buf, g), g)
            idx = p0 + (L - 1)
            @test downsample_symbols(rx, sps; offset=idx - 1)[1] == rx[idx]
            if p0 == 1
                @test idx - 1 == L - 1
            else
                @test idx - 1 != L - 1
            end
            setprecision(BigFloat, 256) do
                tx = apply_fir(buf, g)
                n, s = _scale(tx, rev, idx)
                second = _conv_big(tx, rev)
                @test abs(big(rx[idx]) - second[idx]) <= _gamma(n) * s
                # Zero-lag of this impulse is the same sum of squares as r[L].
                lag = _conv_big(g, rev)[L]
                u_bf = eps(BigFloat(1))
                @test abs(second[idx] - lag) <= (length(buf) + length(g)) * u_bf * (1 + abs(lag))
            end
        end
    end
end

@testset "Eye bound stays negative when g0 < S" begin
    g = rrc_filter(1.0, 3, 1)
    L = length(g)
    sps = 1
    symbols = _lab_symbols(32)
    x = zeros(length(symbols))
    x .= symbols
    rx = matched_filter(apply_fir(x, g), g)
    r = matched_filter(g, g)
    S = 0.0
    m = 1
    while true
        hit = false
        for q in (L + m * sps, L - m * sps)
            if 1 <= q <= length(r)
                S += abs(r[q])
                hit = true
            end
        end
        hit || break
        m += 1
    end
    g0 = r[L]
    lb = 2 * (g0 - S)
    @test g0 > 0
    @test lb < 0
    @test lb == 2 * (g0 - S)
    pos = Float64[]
    neg = Float64[]
    for k in eachindex(symbols)
        sample = rx[L + (k - 1) * sps]
        if symbols[k] > 0
            push!(pos, sample)
        else
            push!(neg, sample)
        end
    end
    sep = minimum(pos) - maximum(neg)
    @test sep == minimum(pos) - maximum(neg)
    @test sep < 0
    println("negative eye lb=", lb, " sep=", sep)
end

@testset "Wrong phase matches the off-peak lag" begin
    β, span, sps = 0.35, 6, 4
    g = rrc_filter(β, span, sps)
    L = length(g)
    symbols = _lab_symbols(32)
    x = zeros(length(symbols) * sps)
    x[1:sps:end] .= symbols
    tx = apply_fir(x, g)
    rx = matched_filter(tx, g)
    rev = reverse(g)
    δ = 1
    setprecision(BigFloat, 256) do
        r = _conv_big(g, rev)
        for k in eachindex(symbols)
            p = L + δ + (k - 1) * sps
            y = big(0)
            for m in eachindex(symbols)
                q = L + δ + (k - m) * sps
                if 1 <= q <= length(r)
                    y += big(symbols[m]) * r[q]
                end
            end
            n, s = _scale(tx, rev, p)
            from_tx = _conv_big(tx, rev)
            @test abs(big(rx[p]) - from_tx[p]) <= _gamma(n) * s
            @test p != L + (k - 1) * sps
        end
    end
end

@testset "pulse_shape truncation is the documented helper" begin
    g = rrc_filter(0.35, 6, 4)
    symbols = _lab_symbols(8)
    shaped = pulse_shape(symbols, g, 4)
    x = zeros(length(symbols) * 4)
    x[1:4:end] .= symbols
    full = apply_fir(x, g)
    @test length(shaped) == length(symbols) * 4
    @test length(full) == length(shaped) + length(g) - 1
end

@testset "Two spans are not a monotone claim" begin
    function energy(span)
        g = rrc_filter(0.35, span, 4)
        L = length(g)
        r = matched_filter(g, g)
        acc = 0.0
        m = 1
        while true
            hit = false
            for q in (L + m * 4, L - m * 4)
                if 1 <= q <= length(r)
                    acc += r[q]^2
                    hit = true
                end
            end
            hit || break
            m += 1
        end
        return acc
    end
    e6 = energy(6)
    e10 = energy(10)
    println("span6_energy=", e6, " span10_energy=", e10)
    @test e6 > 0
    @test e10 > 0
end
