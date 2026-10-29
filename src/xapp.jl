"""
    XApp

Teaching multi-user KPI→action mini xApp: a handful of UEs, each running the
existing single-link teaching loop ([`CommToolkit.Link.simulate_link`](@ref)),
sharing one per-round block budget and one share-weighted mean power budget. A
policy maps cumulative link KPIs (delivered bits, CRC failure rate) onto
discrete actions — an MCS row, a scheduling share, and a power offset — and a
fixed allocation is the baseline it is compared against.

This is a teaching closed loop, not O-RAN: there is no E2 interface, no near-RT
RIC, no fronthaul or functional split, no O-RAN compliance claim, no NR
scheduler, and no NR CQI. The PHY is not rewritten. Every UE is an independent
[`CommToolkit.Link.simulate_link`](@ref) AWGN link; the optional `:rayleigh`
channel is block fading at round granularity with ideal CSI, not a system-level
multi-cell simulator. Multi-user here means independent single links plus the
two shared resource constraints, nothing more.
"""
module XApp

using ..Link: MCS, simulate_link, teaching_punctured_mcs_table,
    teaching_punctured_amc_thresholds, ThresholdAMC, select_mcs
using ..Metrics: snr_db_to_linear

export XAppUE, UEKPI, XAppAction, FixedXAppPolicy, RuleXAppPolicy,
    xapp_actions, observe_kpi!, jain_index, XAppReport, simulate_xapp

# Relative tolerance on the share-weighted mean power check. Shares are
# floating-point weights whose sum can land one ulp above 1; the budget test
# must not fire on that. It is not a way to exceed the budget physically.
const _BUDGET_RTOL = 1e-9

# |h|² of a CN(0,1) draw can underflow the effective Es/N0 to -Inf, which
# simulate_link rejects. A fade below -60 dB is a dropped round for every MCS
# row anyway; P(|h|² < 1e-6) = 1e-6 per draw.
const _FADE_FLOOR_DB = -60.0

@doc raw"""
    XAppUE

One teaching user-equipment descriptor: a name and the base ``E_s/N_0`` of its
own single link, in dB, finite or ``+\infty``, with ``E_s = 1``.

# Fields
- `name`: non-empty label used in reports and tables.
- `esn0_db`: base ``E_s/N_0`` of the UE's AWGN link before any power offset or
  fading draw. Same convention as [`CommToolkit.Link.simulate_link`](@ref).

The descriptor holds no KPI state. Counters live in [`UEKPI`](@ref); the
closed loop is [`simulate_xapp`](@ref). This is not an O-RAN UE model.
"""
struct XAppUE
    name::String
    esn0_db::Float64

    function XAppUE(name::AbstractString, esn0_db::Real)
        !isempty(name) || throw(ArgumentError("UE name must be non-empty"))
        snr = float(esn0_db)
        (snr == Inf || isfinite(snr)) || throw(ArgumentError(
            "Es/N0 must be finite or +Inf, got $esn0_db",
        ))
        return new(String(name), snr)
    end
end

function Base.show(io::IO, ue::XAppUE)
    print(io, "XAppUE(", ue.name, ", Es/N0=", ue.esn0_db, " dB)")
    return nothing
end

@doc raw"""
    XAppAction

The discrete action one policy assigns to one UE for one round.

# Fields
- `mcs`: 1-based row of the MCS table the round is simulated with.
- `share`: fraction of the round's block budget, in ``[0, 1]``. Shares across
  the UEs sum to 1; the runner splits `blocks_per_round` by largest remainder.
- `power_offset_db`: offset added to the UE's base ``E_s/N_0``, in dB, finite.
  ``E_s = 1`` is the 0 dB point; the linear factor ``10^{\mathrm{offset}/10}``
  is what the share-weighted mean power budget charges.

This is a teaching action space. It is not an O-RAN control message and not an
NR scheduling grant.
"""
struct XAppAction
    mcs::Int
    share::Float64
    power_offset_db::Float64

    function XAppAction(mcs::Integer, share::Real, power_offset_db::Real)
        mcs = Int(mcs)
        mcs >= 1 || throw(ArgumentError("mcs index must be ≥ 1, got $mcs"))
        s = float(share)
        (0.0 <= s <= 1.0) || throw(ArgumentError("share must lie in [0, 1], got $share"))
        isnan(s) && throw(ArgumentError("share must not be NaN"))
        o = float(power_offset_db)
        isfinite(o) || throw(ArgumentError(
            "power_offset_db must be finite, got $power_offset_db",
        ))
        return new(mcs, s, o)
    end
end

function Base.show(io::IO, action::XAppAction)
    print(io, "XAppAction(mcs=", action.mcs, ", share=",
          round(action.share; digits=4), ", offset=", action.power_offset_db, " dB)")
    return nothing
end

@doc raw"""
    UEKPI

Cumulative KPI counters of one UE, the only quantities a policy observes.

# Fields
- `name`: copied from [`XAppUE`](@ref).
- `n_blocks`: transport blocks offered so far.
- `n_ack`: blocks whose final decoder word passed the CRC.
- `n_delivered`: blocks with a CRC ACK and a genie-correct payload.
- `n_success_bits`: payload bits in delivered blocks.
- `n_channel_uses`: modulated symbols actually sent, HARQ retransmissions
  included.
- `bler`: `(n_blocks - n_ack) / n_blocks`, and 0 when no block was offered.
- `throughput`: `n_success_bits / n_channel_uses` in bit/use, and 0 when the UE
  sent nothing. At symbol rate 1 and Nyquist bandwidth 1 Hz this is the UE's
  realised spectral efficiency in bit/s/Hz while scheduled.

Counters are cumulative from round 1; there is no sliding window and no
forgetting factor. Built by [`simulate_xapp`](@ref) or by hand for testing a
policy. Not an O-RAN KPI report.
"""
struct UEKPI
    name::String
    n_blocks::Int
    n_ack::Int
    n_delivered::Int
    n_success_bits::Int
    n_channel_uses::Int
    bler::Float64
    throughput::Float64

    function UEKPI(name::AbstractString, n_blocks::Integer, n_ack::Integer,
                   n_delivered::Integer, n_success_bits::Integer,
                   n_channel_uses::Integer)
        !isempty(name) || throw(ArgumentError("UE name must be non-empty"))
        counters = (Int(n_blocks), Int(n_ack), Int(n_delivered),
                    Int(n_success_bits), Int(n_channel_uses))
        all(>=(0), counters) || throw(ArgumentError(
            "KPI counters must be ≥ 0, got $counters",
        ))
        n_delivered <= n_ack <= n_blocks || throw(ArgumentError(
            "need n_delivered ≤ n_ack ≤ n_blocks, got $n_delivered, $n_ack, $n_blocks",
        ))
        bler = n_blocks > 0 ? (n_blocks - n_ack) / n_blocks : 0.0
        throughput = n_channel_uses > 0 ? n_success_bits / n_channel_uses : 0.0
        return new(String(name), counters..., bler, throughput)
    end
end

function Base.show(io::IO, kpi::UEKPI)
    print(io, "UEKPI(", kpi.name, ", blocks=", kpi.n_blocks,
          ", BLER=", round(kpi.bler; sigdigits=3),
          ", T=", round(kpi.throughput; sigdigits=3), " bit/use)")
    return nothing
end

@doc raw"""
    FixedXAppPolicy(; mcs=1, shares=nothing, power_offset_db=0.0)

Baseline allocation: every round assigns the same action to every UE — one
fixed MCS row, a fixed share vector (uniform when `shares` is `nothing`), and
one fixed power offset. It never reads the KPIs. This is the fixed-allocation
baseline the KPI→action [`RuleXAppPolicy`](@ref) is compared against; it is
not an O-RAN policy.

# Fields
- `mcs`: 1-based MCS row, used for every UE and every round.
- `shares`: `nothing` for uniform shares, or a non-negative vector with
  positive sum, normalised to sum 1 at decision time. Its length must match
  the number of UEs.
- `power_offset_db`: constant offset for every UE. The policy does no power
  control, so [`xapp_actions`](@ref) raises `ArgumentError` when the resulting
  share-weighted mean power exceeds `power_budget`.
"""
struct FixedXAppPolicy
    mcs::Int
    shares::Union{Nothing,Vector{Float64}}
    power_offset_db::Float64

    function FixedXAppPolicy(; mcs::Integer=1, shares=nothing,
                             power_offset_db::Real=0.0)
        mcs = Int(mcs)
        mcs >= 1 || throw(ArgumentError("mcs index must be ≥ 1, got $mcs"))
        sh = nothing
        if shares !== nothing
            v = Vector{Float64}(undef, length(shares))
            @inbounds for i in eachindex(shares)
                x = float(shares[i])
                (isfinite(x) && x >= 0) || throw(ArgumentError(
                    "shares must be finite and ≥ 0, got $(shares[i])",
                ))
                v[i] = x
            end
            isempty(v) && throw(ArgumentError("shares must be non-empty"))
            sum(v) > 0 || throw(ArgumentError("shares must have a positive sum"))
            sh = v
        end
        o = float(power_offset_db)
        isfinite(o) || throw(ArgumentError(
            "power_offset_db must be finite, got $power_offset_db",
        ))
        return new(mcs, sh, o)
    end
end

function Base.show(io::IO, policy::FixedXAppPolicy)
    print(io, "FixedXAppPolicy(mcs=", policy.mcs,
          ", shares=", policy.shares === nothing ? "uniform" : policy.shares,
          ", offset=", policy.power_offset_db, " dB)")
    return nothing
end

@doc raw"""
    RuleXAppPolicy(n_ues; target_bler=0.1, offsets_db=[-3.0, 0.0, 3.0],
                   thresholds_db=teaching_punctured_amc_thresholds())

Teaching KPI→action rule table. Every round reads the cumulative [`UEKPI`](@ref)
of every UE and assigns an [`XAppAction`](@ref) by three rules:

1. **Share.** Weight ``w_u = 1/(\texttt{delivered bits}_u + 1)``, normalised.
   A UE that delivered less gets a larger slice of the next round — a
   proportional-fair-like teaching heuristic in the spirit of Tse and
   Viswanath, *Fundamentals of Wireless Communication*, 2005, §6.3, not a 3GPP
   scheduler and not an O-RAN policy. The `+1` guards the first round.
2. **Power offset.** A per-UE notch on the discrete ladder `offsets_db`.
   [`observe_kpi!`](@ref) steps the notch one rung up when the cumulative BLER
   exceeds `target_bler`, one rung down when it falls below
   `target_bler / 10`, and holds otherwise; the notch is clamped to the
   ladder. The initial rung is the offset closest to 0 dB.
3. **MCS.** [`CommToolkit.Link.select_mcs`](@ref) on `thresholds_db` at the
   effective ``E_s/N_0`` = base + current offset. On a static AWGN scenario
   this reads the true base ``E_s/N_0`` (a genie report); under `:rayleigh`
   fading the per-round fade realisation never enters the decision — the
   policy only sees its cumulative KPI consequences, one round late.

The share-weighted mean power ``\sum_u \texttt{share}_u\,10^{\texttt{offset}_u/10}``
must not exceed `power_budget`. The decision projects: while the constraint is
violated, the highest rung in use is lowered one step (ties to the lowest UE
index). A budget that is infeasible even at the minimum rungs raises
`ArgumentError`.

# Fields
- `n_ues`: number of UEs the policy serves; fixed at construction.
- `target_bler`: BLER set point of the offset rule, in ``(0, 1)``.
- `offsets_db`: the power ladder, strictly increasing, finite, non-empty.
- `thresholds_db`: ``E_s/N_0`` switching thresholds, one per MCS row, same
  convention as [`CommToolkit.Link.ThresholdAMC`](@ref). The default is
  [`CommToolkit.Link.teaching_punctured_amc_thresholds`](@ref) and matches the
  default table of [`simulate_xapp`](@ref). They are not NR CQI thresholds.
- `offset_index`: current per-UE rung, mutable state walked by
  [`observe_kpi!`](@ref).

The rule reads cumulative KPIs: early blocks keep contributing to the BLER
forever. That is a teaching simplification, not a sliding-window outer loop.
"""
mutable struct RuleXAppPolicy
    n_ues::Int
    target_bler::Float64
    offsets_db::Vector{Float64}
    thresholds_db::Vector{Float64}
    offset_index::Vector{Int}

    function RuleXAppPolicy(n_ues::Integer;
                            target_bler::Real=0.1,
                            offsets_db::AbstractVector{<:Real}=[-3.0, 0.0, 3.0],
                            thresholds_db::AbstractVector{<:Real}=teaching_punctured_amc_thresholds())
        n_ues = Int(n_ues)
        n_ues >= 1 || throw(ArgumentError("n_ues must be ≥ 1, got $n_ues"))
        (isfinite(target_bler) && 0 < target_bler < 1) || throw(ArgumentError(
            "target_bler must lie in (0, 1), got $target_bler",
        ))
        isempty(offsets_db) && throw(ArgumentError("offsets_db must be non-empty"))
        offsets = Vector{Float64}(undef, length(offsets_db))
        @inbounds for i in eachindex(offsets_db)
            x = float(offsets_db[i])
            isfinite(x) || throw(ArgumentError("offsets_db must be finite, got $x"))
            offsets[i] = x
        end
        for i in 2:length(offsets)
            offsets[i] > offsets[i - 1] || throw(ArgumentError(
                "offsets_db must be strictly increasing, got $(offsets[i - 1]) then $(offsets[i])",
            ))
        end
        thresholds = ThresholdAMC(thresholds_db).thresholds_db
        initial = argmin(abs.(offsets))
        return new(n_ues, float(target_bler), offsets, thresholds,
                   fill(initial, n_ues))
    end
end

function Base.show(io::IO, policy::RuleXAppPolicy)
    print(io, "RuleXAppPolicy(", policy.n_ues, " UEs, target=",
          policy.target_bler, ", offsets=", policy.offsets_db, ")")
    return nothing
end

function _check_budget(power_budget::Real)
    b = float(power_budget)
    (isfinite(b) && b > 0) || throw(ArgumentError(
        "power_budget must be finite and positive, got $power_budget",
    ))
    return b
end

@doc raw"""
    xapp_actions(policy, ues::Vector{XAppUE}, kpis::Vector{UEKPI};
                 power_budget=1.0) -> Vector{XAppAction}

One decision of the KPI→action loop: the actions for the next round, from the
current cumulative KPIs. Pure read of the policy state — the offset ladder of a
[`RuleXAppPolicy`](@ref) moves only in [`observe_kpi!`](@ref).

`power_budget` caps the share-weighted mean linear power
``\sum_u \texttt{share}_u\,10^{\texttt{offset}_u/10}``; the default 1.0 is the
``E_s = 1`` baseline. [`FixedXAppPolicy`](@ref) does no power control: an
offset that breaks the budget raises `ArgumentError`.
[`RuleXAppPolicy`](@ref) lowers the highest rungs in use until the constraint
holds and raises `ArgumentError` only when even the minimum rungs are
infeasible. The check allows a relative ``10^{-9}`` slack for floating-point
share sums.
"""
function xapp_actions(policy::FixedXAppPolicy, ues::AbstractVector{XAppUE},
                      kpis::AbstractVector{UEKPI}; power_budget::Real=1.0)
    n = length(ues)
    length(kpis) == n || throw(ArgumentError(
        "kpis has $(length(kpis)) entries for $n UEs",
    ))
    budget = _check_budget(power_budget)
    w = if policy.shares === nothing
        fill(1.0 / n, n)
    else
        length(policy.shares) == n || throw(ArgumentError(
            "policy has $(length(policy.shares)) shares for $n UEs",
        ))
        policy.shares ./ sum(policy.shares)
    end
    cost = sum(w) * snr_db_to_linear(policy.power_offset_db)
    cost <= budget * (1 + _BUDGET_RTOL) || throw(ArgumentError(
        "FixedXAppPolicy does no power control: offset " *
        "$(policy.power_offset_db) dB puts the share-weighted mean power " *
        "$(cost) above power_budget $budget",
    ))
    return [XAppAction(policy.mcs, w[u], policy.power_offset_db) for u in 1:n]
end

function xapp_actions(policy::RuleXAppPolicy, ues::AbstractVector{XAppUE},
                      kpis::AbstractVector{UEKPI}; power_budget::Real=1.0)
    n = length(ues)
    n == policy.n_ues || throw(ArgumentError(
        "policy was built for $(policy.n_ues) UEs, got $n",
    ))
    length(kpis) == n || throw(ArgumentError(
        "kpis has $(length(kpis)) entries for $n UEs",
    ))
    budget = _check_budget(power_budget)
    # Proportional-fair-like weights: less delivered so far, larger share.
    w = [1.0 / (kpi.n_success_bits + 1) for kpi in kpis]
    shares = w ./ sum(w)
    offsets = policy.offsets_db
    notch = copy(policy.offset_index)
    function mean_power()
        c = 0.0
        @inbounds for u in 1:n
            c += shares[u] * snr_db_to_linear(offsets[notch[u]])
        end
        return c
    end
    while mean_power() > budget * (1 + _BUDGET_RTOL)
        movable = findall(>(1), notch)
        isempty(movable) && throw(ArgumentError(
            "power_budget $budget is infeasible even at the minimum rungs " *
            "of offsets_db $offsets; raise the budget or lower the ladder",
        ))
        # Lower the highest rung in use; ties go to the lowest UE index.
        u = movable[argmax(offsets[notch[movable]])]
        notch[u] -= 1
    end
    amc = ThresholdAMC(policy.thresholds_db)
    return [XAppAction(select_mcs(amc, ues[u].esn0_db + offsets[notch[u]]),
                       shares[u], offsets[notch[u]]) for u in 1:n]
end

@doc raw"""
    observe_kpi!(policy::RuleXAppPolicy, ue::Integer, kpi::UEKPI) -> Int
    observe_kpi!(policy::FixedXAppPolicy, ue::Integer, kpi::UEKPI) -> nothing

Feed one UE's cumulative KPI back to the rule table. The `RuleXAppPolicy`
method walks that UE's power-offset rung: one rung up when `kpi.bler` exceeds
`policy.target_bler`, one rung down when it is below `target_bler / 10`, hold
otherwise, clamped to `offsets_db`. A KPI with `n_blocks == 0` is no
observation and does not move the rung. The return is the current rung index.
The `FixedXAppPolicy` method is a no-op, so a fixed baseline never reacts to
KPIs.

[`simulate_xapp`](@ref) calls this once per UE per round, after that round's
blocks are counted, with the cumulative KPI.
"""
function observe_kpi!(policy::RuleXAppPolicy, ue::Integer, kpi::UEKPI)
    u = Int(ue)
    1 <= u <= policy.n_ues || throw(ArgumentError(
        "UE index $u outside 1:$(policy.n_ues)",
    ))
    if kpi.n_blocks >= 1
        if kpi.bler > policy.target_bler
            policy.offset_index[u] = min(policy.offset_index[u] + 1,
                                         length(policy.offsets_db))
        elseif kpi.bler < policy.target_bler / 10
            policy.offset_index[u] = max(policy.offset_index[u] - 1, 1)
        end
    end
    return policy.offset_index[u]
end

observe_kpi!(policy::FixedXAppPolicy, ue::Integer, kpi::UEKPI) = nothing

@doc raw"""
    jain_index(x::AbstractVector{<:Real}) -> Float64

Jain's fairness index

```math
J(x) = \frac{\bigl(\sum_u x_u\bigr)^2}{n \sum_u x_u^2}
```

of a non-negative allocation (Jain, Chiu, and Hawe, DEC-TR-301, 1984). Equal
shares give 1; one UE holding everything gives ``1/n``. An all-zero vector
returns `NaN` — fairness is undefined when nothing was delivered. Negative or
non-finite entries raise `ArgumentError`.

[`simulate_xapp`](@ref) reports it over the per-UE delivered payload bits.
"""
function jain_index(x::AbstractVector{<:Real})
    n = length(x)
    n >= 1 || throw(ArgumentError("jain_index needs a non-empty vector"))
    s1 = 0.0
    s2 = 0.0
    @inbounds for i in eachindex(x)
        v = float(x[i])
        (isfinite(v) && v >= 0) || throw(ArgumentError(
            "allocations must be finite and ≥ 0, got $(x[i])",
        ))
        s1 += v
        s2 += v * v
    end
    s1 == 0 && return NaN
    return s1 * s1 / (n * s2)
end

@doc raw"""
    XAppReport

Aggregate of one [`simulate_xapp`](@ref) run.

# Fields
- `ues`: the UE descriptors, in order.
- `kpis`: cumulative per-UE [`UEKPI`](@ref) at the end of the run.
- `n_rounds`, `blocks_per_round`: the round count and the per-round block
  budget that were simulated.
- `total_blocks`, `total_ack`, `total_delivered`: sums over UEs and rounds.
- `total_success_bits`, `total_channel_uses`: delivered payload bits and
  symbols sent, summed over UEs and rounds.
- `system_bler`: `1 - total_ack / total_blocks`. 0 when nothing was offered.
- `system_throughput`: `total_success_bits / total_channel_uses` in bit/use.
  0 when nothing was sent. At symbol rate 1 and Nyquist bandwidth 1 Hz this is
  the aggregate spectral efficiency of the scheduled links in bit/s/Hz.
- `jain`: [`jain_index`](@ref) of the per-UE delivered payload bits. `NaN`
  when no UE delivered anything.
- `last_actions`: the actions of the final round, one per UE.

Not an O-RAN report format.
"""
struct XAppReport
    ues::Vector{XAppUE}
    kpis::Vector{UEKPI}
    n_rounds::Int
    blocks_per_round::Int
    total_blocks::Int
    total_ack::Int
    total_delivered::Int
    total_success_bits::Int
    total_channel_uses::Int
    system_bler::Float64
    system_throughput::Float64
    jain::Float64
    last_actions::Vector{XAppAction}
end

function Base.show(io::IO, report::XAppReport)
    print(io, "XAppReport(", length(report.ues), " UEs × ", report.n_rounds,
          " rounds, T=", round(report.system_throughput; sigdigits=3),
          " bit/use, BLER=", round(report.system_bler; sigdigits=3),
          ", Jain=", round(report.jain; sigdigits=3), ")")
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", report::XAppReport)
    println(io, "XAppReport")
    println(io, "  rounds: ", report.n_rounds, " × ", report.blocks_per_round, " blocks")
    println(io, "  per-UE:  name    blocks  delivered    BLER   T (bit/use)")
    for kpi in report.kpis
        println(io, "           ", rpad(kpi.name, 8), lpad(kpi.n_blocks, 6),
                lpad(kpi.n_delivered, 10), "  ", round(kpi.bler; sigdigits=3),
                "   ", round(kpi.throughput; sigdigits=3))
    end
    println(io, "  system throughput (bit/use): ", report.system_throughput)
    println(io, "  system BLER: ", report.system_bler)
    println(io, "  Jain fairness (delivered bits): ", report.jain)
    print(io, "  last actions: ", report.last_actions)
    return nothing
end

"""Largest-remainder split of `total` blocks by `shares`; ties to the lower index."""
function _split_blocks(total::Int, shares::Vector{Float64})
    raw = total .* shares
    counts = floor.(Int, raw)
    leftover = total - sum(counts)
    frac = raw .- counts
    for _ in 1:leftover
        u = argmax(frac)
        counts[u] += 1
        frac[u] = -1.0
    end
    return counts
end

@doc raw"""
    simulate_xapp(ues::Vector{XAppUE}, policy, n_rounds::Integer;
                  kwargs...) -> XAppReport

Teaching multi-user KPI→action closed loop. Every round: the policy reads the
cumulative [`UEKPI`](@ref)s and assigns an [`XAppAction`](@ref) per UE
([`xapp_actions`](@ref)); the round's `blocks_per_round` transport blocks are
split by the shares (largest remainder, ties to the lower UE index); each UE
with at least one block runs the existing single-link
[`CommToolkit.Link.simulate_link`](@ref) at

```math
E_s/N_0 = \texttt{base}_u + \texttt{offset}_u \; (+\; 10\log_{10}|h_u|^2
\ \text{when}\ \texttt{channel} = \texttt{:rayleigh})
```

with the action's MCS row; the counters accumulate; the rule table observes
the new cumulative KPIs ([`observe_kpi!`](@ref)).

Multi-user here means independent single links under two shared constraints:
the per-round block budget and the share-weighted mean power budget
``\sum_u \texttt{share}_u\,10^{\texttt{offset}_u/10} \le \texttt{power\_budget}``
(default 1.0, the ``E_s = 1`` baseline). There is no inter-UE interference
model, no cell geometry, and no HARQ timeline across UEs.

# Arguments
- `ues`: UE descriptors, at least one.
- `policy`: [`FixedXAppPolicy`](@ref) (the fixed-allocation baseline) or
  [`RuleXAppPolicy`](@ref) (the KPI→action rule table).
- `n_rounds`: number of scheduling rounds.

# Keywords
- `blocks_per_round`: transport-block budget of one round, shared by the UEs.
  Default 30.
- `table`: MCS list, default
  [`CommToolkit.Link.teaching_punctured_mcs_table`](@ref) (all convolutional,
  so `harq = :ir` works on every row). A `RuleXAppPolicy`'s `thresholds_db`
  must match its length.
- `max_tx`, `harq`, `approx`: forwarded to
  [`CommToolkit.Link.simulate_link`](@ref) unchanged. `harq = :ir` keeps that
  function's checks (convolutional rows, complement length). The receiver is
  always the default `:viterbi`; the BICM receiver stays on the single-link
  entry point.
- `channel`: `:awgn` (default) or `:rayleigh`. `:rayleigh` draws one
  ``h_u \sim \mathcal{CN}(0,1)`` per UE per round and adds
  ``10\log_{10}|h_u|^2`` (clamped at −60 dB) to that round's ``E_s/N_0``.
  With ideal CSI the faded link conditioned on ``h`` *is* the AWGN link at
  ``|h|^2 E_s/N_0`` — the same conditional equivalence as
  `rayleigh_channel(..., true)` — so this is block fading at round
  granularity, not a new PHY. The fade realisation never enters the policy;
  it reacts one round late through the KPIs alone.
- `power_budget`: share-weighted mean power cap, finite and positive.

# Random stream
Global RNG, like [`CommToolkit.Link.simulate_link`](@ref). With
`channel = :rayleigh` each round first draws one ``h`` per UE in index order
(whether or not the UE gets blocks), then the per-UE `simulate_link` calls
follow in index order. With `:awgn` no fading draws happen. Seed at the call
site (`Random.seed!`) to reproduce a run; two policies compared on the same
seed see the same stream layout.

This is a teaching mini xApp, not O-RAN: no E2 interface, no near-RT RIC, no
fronthaul, no NR scheduler, no NR CQI, and no O-RAN compliance claim. The MCS
thresholds are the teaching ones; the fairness heuristic is not a 3GPP
scheduler.

# References
Jain, Chiu, and Hawe, DEC-TR-301, 1984 (fairness index). Tse and Viswanath,
*Fundamentals of Wireless Communication*, 2005, §6.3 (proportional fairness;
the share rule is only a like-flavoured teaching heuristic). Switching
thresholds: Goldsmith, 2005, §9.3, via
[`CommToolkit.Link.teaching_punctured_amc_thresholds`](@ref). Link chain and
HARQ: [`CommToolkit.Link.simulate_link`](@ref).
"""
function simulate_xapp(ues::AbstractVector{XAppUE},
                       policy::Union{FixedXAppPolicy,RuleXAppPolicy},
                       n_rounds::Integer;
                       blocks_per_round::Integer=30,
                       table::AbstractVector{MCS}=teaching_punctured_mcs_table(),
                       max_tx::Integer=1,
                       harq::Symbol=:chase,
                       approx::Symbol=:exact,
                       channel::Symbol=:awgn,
                       power_budget::Real=1.0)
    n_ues = length(ues)
    n_ues >= 1 || throw(ArgumentError("ues must be non-empty"))
    n_rounds = Int(n_rounds)
    n_rounds >= 1 || throw(ArgumentError("n_rounds must be ≥ 1, got $n_rounds"))
    blocks_per_round = Int(blocks_per_round)
    blocks_per_round >= 1 || throw(ArgumentError(
        "blocks_per_round must be ≥ 1, got $blocks_per_round",
    ))
    isempty(table) && throw(ArgumentError("MCS table must be non-empty"))
    max_tx = Int(max_tx)
    max_tx >= 1 || throw(ArgumentError("max_tx must be ≥ 1, got $max_tx"))
    harq === :chase || harq === :ir || throw(ArgumentError(
        "harq must be :chase or :ir, got $(repr(harq))",
    ))
    approx === :exact || approx === :maxlog || throw(ArgumentError(
        "approx must be :exact or :maxlog, got $(repr(approx))",
    ))
    channel === :awgn || channel === :rayleigh || throw(ArgumentError(
        "channel must be :awgn or :rayleigh, got $(repr(channel))",
    ))
    budget = _check_budget(power_budget)
    if policy isa RuleXAppPolicy
        policy.n_ues == n_ues || throw(ArgumentError(
            "policy was built for $(policy.n_ues) UEs, got $n_ues",
        ))
        length(policy.thresholds_db) == length(table) || throw(ArgumentError(
            "policy has $(length(policy.thresholds_db)) thresholds for a " *
            "table of $(length(table))",
        ))
    else
        1 <= policy.mcs <= length(table) || throw(ArgumentError(
            "fixed mcs index $(policy.mcs) outside 1:$(length(table))",
        ))
        (policy.shares === nothing || length(policy.shares) == n_ues) ||
            throw(ArgumentError(
                "policy has $(length(policy.shares)) shares for $n_ues UEs",
            ))
    end

    n_blocks = zeros(Int, n_ues)
    n_ack = zeros(Int, n_ues)
    n_delivered = zeros(Int, n_ues)
    n_success_bits = zeros(Int, n_ues)
    n_channel_uses = zeros(Int, n_ues)
    actions = XAppAction[]

    for _round in 1:n_rounds
        kpis = UEKPI[UEKPI(ues[u].name, n_blocks[u], n_ack[u], n_delivered[u],
                           n_success_bits[u], n_channel_uses[u]) for u in 1:n_ues]
        actions = xapp_actions(policy, ues, kpis; power_budget=budget)
        fades_db = zeros(n_ues)
        if channel === :rayleigh
            # One CN(0,1) draw per UE per round, index order, before any
            # per-UE link call, whether or not the UE gets blocks.
            for u in 1:n_ues
                h = (randn() + im * randn()) / sqrt(2)
                fades_db[u] = max(10 * log10(abs2(h)), _FADE_FLOOR_DB)
            end
        end
        counts = _split_blocks(blocks_per_round, [a.share for a in actions])
        for u in 1:n_ues
            counts[u] >= 1 || continue
            eff_db = ues[u].esn0_db + actions[u].power_offset_db + fades_db[u]
            report = simulate_link(eff_db, counts[u]; table=table,
                                   mcs=actions[u].mcs, max_tx=max_tx,
                                   harq=harq, approx=approx)
            n_blocks[u] += report.n_blocks
            n_ack[u] += report.n_ack
            n_delivered[u] += report.n_delivered
            n_success_bits[u] += report.n_success_bits
            n_channel_uses[u] += report.n_channel_uses
        end
        if policy isa RuleXAppPolicy
            newkpis = UEKPI[UEKPI(ues[u].name, n_blocks[u], n_ack[u],
                                  n_delivered[u], n_success_bits[u],
                                  n_channel_uses[u]) for u in 1:n_ues]
            for u in 1:n_ues
                observe_kpi!(policy, u, newkpis[u])
            end
        end
    end

    kpis = UEKPI[UEKPI(ues[u].name, n_blocks[u], n_ack[u], n_delivered[u],
                       n_success_bits[u], n_channel_uses[u]) for u in 1:n_ues]
    total_blocks = sum(n_blocks)
    total_ack = sum(n_ack)
    total_bits = sum(n_success_bits)
    total_uses = sum(n_channel_uses)
    return XAppReport(
        collect(ues), kpis, n_rounds, blocks_per_round,
        total_blocks, total_ack, sum(n_delivered), total_bits, total_uses,
        total_blocks > 0 ? 1 - total_ack / total_blocks : 0.0,
        total_uses > 0 ? total_bits / total_uses : 0.0,
        jain_index(n_success_bits),
        actions,
    )
end

end # module XApp
