module Viewer

using GLMakie
using Zarr

const U = parentmodule(@__MODULE__)

const MAP_BACK = "Pr(map-back)"
const EX_RATE = "ex-map rate | o.l."
const LOG_RATIO = "log2(map-back / ex-map rate)"
const METRICS = [MAP_BACK, EX_RATE, LOG_RATIO]
const PREV, NEXT = "\u25C0", "\u25B6"

function open_store(path::AbstractString; read_lengths=nothing)
    isdir(path) || error("Zarr store not found: $path")
    g = zopen(path; consolidated=isfile(joinpath(path, ".zmetadata")))
    haskey(g.arrays, "n_correct") || error("$path has no count arrays; re-run `unimap score` with this version")
    Ls = Int.(g.attrs["read_lengths"])
    iLs = read_lengths === nothing ? collect(eachindex(Ls)) :
          [something(findfirst(==(L), Ls), 0) for L in read_lengths]
    any(==(0), iLs) && error("read lengths $read_lengths not all in store ($Ls)")
    loci = g.attrs["loci"]
    return (; g, loci, ids=String[l["locus_id"] for l in loci], Ls, iLs,
            Fs=Int.(g.attrs["fragment_lengths"]), cuts=Int.(g.attrs["mapq_cutoffs"]),
            pseudocount=Float64(get(g.attrs, "log_ratio_pseudocount", U.LOG_RATIO_PSEUDOCOUNT)))
end

function locus_index(st, locus)
    locus === nothing && return 1
    k = findfirst(==(locus), st.ids)
    k === nothing && error("locus '$locus' not in store")
    return k
end

nlen(st, k) = st.loci[k]["end0"] - st.loci[k]["start0"]

"""Counts for one locus, restricted to the shown read lengths: ns[x,F,L], nc/ne[x,F,cutoff,L]."""
function load_block(st, k)
    rows = st.loci[k]["offset"]+1:st.loci[k]["offset"]+nlen(st, k)
    return (; ns=st.g["n_simulated"][rows, :, :][:, :, st.iLs],
            nc=st.g["n_correct"][rows, :, :, :][:, :, :, st.iLs],
            ne=st.g["n_ectopic"][rows, :, :, :][:, :, :, st.iLs], pseudocount=st.pseudocount)
end

function locus_label(st, k)
    l = st.loci[k]
    return "$(l["chrom"]):$(l["start0"] + 1)-$(l["end0"]) ($(l["strand"]))    $(nlen(st, k)) bp    locus $k of $(length(st.ids))"
end

function metric_array(m, ns, nc, ne, pc)
    m == MAP_BACK && return U.map_back.(nc, ns)
    m == EX_RATE && return U.ectopic_rate.(ne, ns)
    return U.log2_ratio.(nc, ne; pseudocount=pc)
end

"""Metric at base `i` for read-length column `col`: matrix over (F, cutoff)."""
function metric_at_base(b, i, col, m)
    i = clamp(i, 1, size(b.ns, 1))
    return metric_array(m, b.ns[i, :, col], b.nc[i, :, :, col], b.ne[i, :, :, col], b.pseudocount)
end

"""Metric along the locus at cutoff index `ic` for read-length column `col`: matrix over (position, F)."""
metric_at_cutoff(b, ic, col, m) =
    metric_array(m, b.ns[:, :, col], b.nc[:, :, ic, col], b.ne[:, :, ic, col], b.pseudocount)

finite_max(f, xs) = maximum((v for v in (f(x) for x in xs) if isfinite(v)); init=0f0)

"""z range and colormap for a metric, fixed over the whole locus so values are comparable along it."""
function metric_scale(b, m)
    m == MAP_BACK && return (0f0, 1f0), :viridis
    ns = reshape(b.ns, size(b.ns, 1), size(b.ns, 2), 1, size(b.ns, 3))
    if m == EX_RATE
        return (0f0, max(1f0, finite_max(identity, U.ectopic_rate.(b.ne, ns)))), :viridis
    end
    M = max(1f0, finite_max(abs, U.log2_ratio.(b.nc, b.ne; pseudocount=b.pseudocount)))
    return (-M, M), :RdBu
end

"""Locus dropdown with prev/next buttons, metric dropdown and optionally a MAPQ cutoff dropdown."""
function controls!(fig, st, k0, ncols; cutoff=nothing)
    ctrl = fig[1, 1:ncols] = GridLayout(tellwidth=false)
    prev = Button(ctrl[1, 1], label=PREV)
    menu = Menu(ctrl[1, 2], options=st.ids, default=st.ids[k0], width=260)
    next = Button(ctrl[1, 3], label=NEXT)
    metric_menu = Menu(ctrl[1, 4], options=METRICS, default=METRICS[1], width=230)
    col, cutoff_menu = 5, nothing
    if cutoff !== nothing
        cutoff in st.cuts || error("MAPQ cutoff $cutoff is not in the store (available: $(join(st.cuts, ", ")))")
        Label(ctrl[1, 5], "MAPQ >=")
        cutoff_menu = Menu(ctrl[1, 6], options=string.(st.cuts), default=string(cutoff), width=70)
        col = 7
    end
    k = Observable(k0)
    Label(ctrl[1, col], lift(kk -> locus_label(st, kk), k), halign=:left, tellwidth=false)
    on(menu.selection) do id
        kk = findfirst(==(id), st.ids)
        kk === nothing || kk == k[] || (k[] = kk)
    end
    on(_ -> (menu.i_selected[] = mod1(k[] - 1, length(st.ids))), prev.clicks)
    on(_ -> (menu.i_selected[] = mod1(k[] + 1, length(st.ids))), next.clicks)
    return (; k, menu, metric_menu, cutoff_menu)
end

cutoff_index(st, c) = findfirst(==(parse(Int, c.cutoff_menu.selection[])), st.cuts)

# GLMakie uploads a surface as one texture; Apple's OpenGL caps texture sides at 16384.
const MAX_SURFACE_POINTS = 8192

worst(m) = m == EX_RATE ? max : min
worst_word(m) = m == EX_RATE ? "max" : "min"

"""
    bin_positions(xs, vals, m) -> (bin centres, binned values, bin width)

If there are more than `MAX_SURFACE_POINTS` positions, reduce `vals` (position x F) to
equal-width bins keeping the worst finite value of metric `m` in each bin (lowest
map-back or log ratio, highest ex-map rate), so dips are not averaged away.
"""
function bin_positions(xs, vals::AbstractMatrix, m)
    n = size(vals, 1)
    w = cld(n, MAX_SURFACE_POINTS)
    w == 1 && return collect(xs), vals, 1
    f, nb = worst(m), cld(n, w)
    out = Matrix{Float32}(undef, nb, size(vals, 2))
    centres = Vector{Float64}(undef, nb)
    for bi in 1:nb
        r = (bi-1)*w+1:min(bi * w, n)
        centres[bi] = (xs[first(r)] + xs[last(r)]) / 2
        for j in axes(vals, 2)
            v = filter(!isnan, view(vals, r, j))
            out[bi, j] = isempty(v) ? NaN32 : reduce(f, v)
        end
    end
    return centres, out, w
end

"""
    build_viewer(zarr_path; read_lengths=nothing, locus=nothing)

Per-base view. One 3-D axis per read length shows, at a single base `s`, the chosen metric
as a surface over (fragment length x MAPQ cutoff), linearly interpolated between grid
points. The slider (or Left/Right for +/-1 bp, Up/Down for +/-100 bp) moves `s`.
"""
function build_viewer(path::AbstractString; read_lengths=nothing, locus=nothing)
    st = open_store(path; read_lengths)
    k0, nL = locus_index(st, locus), length(st.iLs)
    fig = Figure(size=(480 * nL + 140, 700))
    c = controls!(fig, st, k0, nL)
    metric = c.metric_menu.selection
    block = Observable(load_block(st, k0))
    slider = Slider(fig[3, 1:nL], range=1:nlen(st, k0), startvalue=1)
    s = lift((b, v) -> clamp(v, 1, size(b.ns, 1)), block, slider.value)
    scale = lift(metric_scale, block, metric)
    zlims, cmap = lift(first, scale), lift(last, scale)
    Fs, cuts = st.Fs, st.cuts

    for (col, t) in enumerate(st.iLs)
        ax = Axis3(fig[2, col]; title="L = $(st.Ls[t])", xlabel="fragment length (bp)",
                   ylabel="MAPQ cutoff", zlabel=metric[], azimuth=-0.35pi, elevation=0.12pi)
        on(m -> (ax.zlabel = m), metric)
        z = lift((b, i, m) -> metric_at_base(b, i, col, m), block, s, metric)
        surface!(ax, Fs, cuts, z; colormap=cmap, colorrange=zlims)
        wireframe!(ax, Fs, cuts, z; color=(:black, 0.6), linewidth=0.6)
        scatter!(ax, lift(zz -> vec([Point3f(F, q, zz[i, j]) for (i, F) in enumerate(Fs), (j, q) in enumerate(cuts)]), z);
                 color=:black, markersize=5)
        on(zlims) do (lo, hi)
            limits!(ax, minimum(Fs), maximum(Fs), minimum(cuts), maximum(cuts), lo, hi)
        end
        notify(zlims)
    end
    Colorbar(fig[2, nL+1], colormap=cmap, limits=zlims)
    colgap!(fig.layout, nL, 40)

    Label(fig[4, 1:nL], lift(c.k, s) do kk, i
        l = st.loci[kk]
        "s = $i / $(nlen(st, kk))    $(l["chrom"]):$(l["start0"] + i) (1-based)    Left/Right +/-1 bp, Up/Down +/-100 bp"
    end)

    on(c.k) do kk
        block[] = load_block(st, kk)
        slider.range[] = 1:nlen(st, kk)
        set_close_to!(slider, clamp(slider.value[], 1, nlen(st, kk)))
    end
    on(events(fig).keyboardbutton) do ev
        ev.action in (Keyboard.press, Keyboard.repeat) || return
        d = ev.key == Keyboard.right ? 1 : ev.key == Keyboard.left ? -1 :
            ev.key == Keyboard.up ? 100 : ev.key == Keyboard.down ? -100 : 0
        d == 0 || set_close_to!(slider, clamp(slider.value[] + d, 1, size(block[].ns, 1)))
    end

    return (; figure=fig, menu=c.menu, metric_menu=c.metric_menu, slider, block)
end

"""
    build_track_viewer(zarr_path; read_lengths=nothing, locus=nothing, mapq_cutoff=2)

Along-locus view. One 3-D axis per read length shows the chosen metric at a fixed MAPQ
cutoff (changeable from a dropdown) as a surface over (position along the locus x fragment
length). Undefined values (NaN) are left as holes. Loci longer than `MAX_SURFACE_POINTS`
are shown in bins holding the worst value per bin (see `bin_positions`).
"""
function build_track_viewer(path::AbstractString; read_lengths=nothing, locus=nothing, mapq_cutoff::Integer=2)
    st = open_store(path; read_lengths)
    k0, nL = locus_index(st, locus), length(st.iLs)
    fig = Figure(size=(560 * nL + 140, 720))
    c = controls!(fig, st, k0, nL; cutoff=mapq_cutoff)
    axs = [Axis3(fig[2, col]; title="L = $(st.Ls[t])", xlabel="position (bp, 1-based)",
                 ylabel="fragment length (bp)", zlabel=METRICS[1], azimuth=-0.3pi, elevation=0.15pi)
           for (col, t) in enumerate(st.iLs)]
    cmap, lims = Observable{Any}(:viridis), Observable((0f0, 1f0))
    Colorbar(fig[2, nL+1], colormap=cmap, limits=lims)
    colgap!(fig.layout, nL, 40)
    block = Ref(load_block(st, k0))

    function redraw()
        b, m, ic = block[], c.metric_menu.selection[], cutoff_index(st, c)
        (lo, hi), cm = metric_scale(b, m)
        l = st.loci[c.k[]]
        for (col, ax) in enumerate(axs)
            empty!(ax)
            xs, vals, w = bin_positions(l["start0"]+1:l["end0"], metric_at_cutoff(b, ic, col, m), m)
            surface!(ax, xs, st.Fs, vals; colormap=cm, colorrange=(lo, hi))
            ax.zlabel = m
            ax.title = "L = $(st.Ls[st.iLs[col]])" * (w > 1 ? "  ($w-bp bins, $(worst_word(m)) per bin)" : "")
            limits!(ax, first(xs), last(xs), minimum(st.Fs), maximum(st.Fs), lo, hi)
        end
        cmap[], lims[] = cm, (lo, hi)
    end

    on(kk -> (block[] = load_block(st, kk); redraw()), c.k)
    on(_ -> redraw(), c.metric_menu.selection)
    on(_ -> redraw(), c.cutoff_menu.selection)
    redraw()
    return (; figure=fig, menu=c.menu, metric_menu=c.metric_menu, cutoff_menu=c.cutoff_menu, block)
end

"""
    build_boxplot_viewer(zarr_path; read_lengths=nothing, locus=nothing, mapq_cutoff=2)

Distribution view. One panel per read length with a box plot, per fragment length, of the
chosen metric across all bases of the locus at a fixed MAPQ cutoff (changeable from a
dropdown). Undefined values (NaN) are excluded and counted in the panel title.
"""
function build_boxplot_viewer(path::AbstractString; read_lengths=nothing, locus=nothing, mapq_cutoff::Integer=2)
    st = open_store(path; read_lengths)
    k0, nL, nF = locus_index(st, locus), length(st.iLs), length(st.Fs)
    fig = Figure(size=(460 * nL + 80, 620))
    c = controls!(fig, st, k0, nL; cutoff=mapq_cutoff)
    axs = [Axis(fig[2, col]; xlabel="fragment length (bp)", ylabel=METRICS[1],
                xticks=(1:nF, string.(st.Fs)), xticklabelrotation=pi / 4)
           for col in 1:nL]
    linkyaxes!(axs...)
    block = Ref(load_block(st, k0))

    function redraw()
        b, m, ic = block[], c.metric_menu.selection[], cutoff_index(st, c)
        (lo, hi), _ = metric_scale(b, m)
        for (col, ax) in enumerate(axs)
            empty!(ax)
            vals = metric_at_cutoff(b, ic, col, m)
            xs, ys = repeat(1:nF, inner=size(vals, 1)), vec(vals)
            ok = isfinite.(ys)
            any(ok) && boxplot!(ax, xs[ok], ys[ok]; color=(:steelblue, 0.7), markersize=3,
                                whiskerwidth=0.5, mediancolor=:black)
            m == LOG_RATIO && hlines!(ax, [0]; color=:gray, linestyle=:dash)
            nundef = count(!, ok)
            ax.title = "L = $(st.Ls[st.iLs[col]])" * (nundef > 0 ? "  ($(round(100nundef / length(ok); digits=1))% undefined)" : "")
            ax.ylabel = m
            pad = 0.04f0 * (hi - lo)
            ylims!(ax, lo - pad, hi + pad)
            xlims!(ax, 0.4, nF + 0.6)
        end
    end

    on(kk -> (block[] = load_block(st, kk); redraw()), c.k)
    on(_ -> redraw(), c.metric_menu.selection)
    on(_ -> redraw(), c.cutoff_menu.selection)
    redraw()
    return (; figure=fig, menu=c.menu, metric_menu=c.metric_menu, cutoff_menu=c.cutoff_menu, block)
end

const BUILDERS = Dict(:base => build_viewer, :track => build_track_viewer, :boxplot => build_boxplot_viewer)

function build(path::AbstractString; plot::Symbol=:base, mapq_cutoff::Integer=2, kw...)
    haskey(BUILDERS, plot) || error("plot must be one of: base, track, boxplot (got $plot)")
    return plot == :base ? build_viewer(path; kw...) : BUILDERS[plot](path; mapq_cutoff, kw...)
end

function open_viewer(path::AbstractString; kw...)
    v = build(path; kw...)
    display(v.figure)
    return v
end

"""Open the viewer in a window and block until it is closed."""
function run_viewer(path::AbstractString; kw...)
    v = build(path; kw...)
    wait(display(v.figure))
    return nothing
end

end
