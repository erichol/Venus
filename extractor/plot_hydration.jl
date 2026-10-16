#!/usr/bin/env julia

using MAT
using Plots

gr()  # 3D backend

# ==============================
# Input .mat file
# ==============================
matfile = length(ARGS) >= 1 ? ARGS[1] : "database_VBasalt.mat"
println("Reading MAT file: $matfile")

vars = matread(matfile)

rhosolid  = vars["rhosolid"]      # (np × nt × nwt)
rhosys    = vars["rhosys"]        # (np × nt × nwt)
tkj       = vec(vars["tkj"])      # T axis (K), length nt
pbi       = vec(vars["pbi"])      # P axis (GPa), length np
wth20k    = vec(vars["wth20k"])   # hydration axis, length nwt

np, nt, nwt = size(rhosolid)
println("Grid sizes: np = $np (P), nt = $nt (T), nwt = $nwt (hydrations)")

# ==============================
# Threshold for "no free water"
# ==============================
ϵ = 0.1    # kg/m³  <-- lower threshold
println("Using tolerance ϵ = $ϵ kg/m³ for 'no free water'.")

# ==============================
# Collect all free-water points
# ==============================
allP = Float64[]
allT = Float64[]
allZ = Float64[]   # hydration coordinate (H₂O_mol)

for iw in 1:nwt
    hval = wth20k[iw]

    ρs   = rhosolid[:, :, iw]   # (np, nt)
    ρsys = rhosys[:, :, iw]
    Δρ   = ρs .- ρsys

    if all(isnan, Δρ)
        @warn "All NaN for hydration = $hval, skipping"
        continue
    end

    # mask: true => free water (densities differ enough)
    freewater = abs.(Δρ) .>= ϵ      # (np, nt)

    # Build P,T grids to match (np, nt), then transpose to (nt, np)
    Pmat = repeat(pbi', nt, 1)      # (nt, np)
    Tmat = repeat(tkj, 1, np)       # (nt, np)
    fwT  = freewater'               # (nt, np)

    # Flatten and keep only free-water points
    idx = vec(fwT .== true)
    if !any(idx)
        @warn "No free water points for hydration = $hval at ϵ = $ϵ"
        continue
    end

    append!(allP, vec(Pmat)[idx])
    append!(allT, vec(Tmat)[idx])
    append!(allZ, fill(hval, count(idx)))
end

if isempty(allP)
    error("No free water points found for any hydration at ϵ = $ϵ")
end

# ==============================
# 3D scatter: points = free water
# color = hydration level
# ==============================
plt = scatter(
    allP, allT, allZ,
    xlab   = "P (GPa)",
    ylab   = "T (K)",
    zlab   = "H₂O_mol (hydration)",
    title  = "Free water presence across P–T–hydration",
    zcolor = allZ,                     # color by hydration
    colorbar_title = "H₂O_mol",
    m      = (:circle, 3, 0.5),        # (shape, size, alpha)
    markerstrokewidth = 0,
    camera = (50, 40),
    grid   = false,
    background_color = :white,
)

xlims!(minimum(pbi), maximum(pbi))
ylims!(minimum(tkj), maximum(tkj))
zlims!(minimum(wth20k), maximum(wth20k))

# Static PNG
outpng = "PT_freewater_3D_scatter.png"
println("Saving 3D free-water scatter plot to $outpng")
savefig(plt, outpng)

# ==============================
# Optional: spinning GIF
# ==============================
anim = @animate for ang in 0:5:355
    plot(plt, camera=(ang, 30))
end

outgif = "PT_freewater_3D_scatter_spin.gif"
println("Saving spinning GIF to $outgif")
gif(anim, outgif, fps=20)

println("Done.")
