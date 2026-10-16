using MAGEMin_C, DelimitedFiles, GLMakie, Statistics, ForwardDiff, QuadGK
const R = 8314510 # Gas constant in Pa cm³ K⁻¹ mol⁻¹
const MM = [60.0843, 40.304, 71.844, 159.688, 56.077, 101.961, 79.865, 70.937, 94.196, 61.979, 18.015]

@inline function wt2mol(X)
    ans = copy(X); ans ./= MM; ans ./= sum(ans); ans .*= 1e2
    return ans
end

function Hwt2mol(Hwt, X)
    Hmol = zeros(length(Hwt))
    Xu = copy(X)
    for i in eachindex(Hwt)
        Xu[end] = Hwt[i]; 
        Hmol[i] = wt2mol(Xu)[end]
    end
    return Hmol
end

function XH2Osolid_sH(out)
    XH2Osolid = 0.0
    H = 0.0; solid_n = 0.0
    # Solution Phases
    for entry in 1:length(out.SS_vec)
        out.ph[entry] == "fl" && continue
        XH2Osolid += out.ph_frac[entry] * out.SS_vec[entry].Comp[6] # Xmol_phase * H₂O content
        solid_n += out.ph_frac[entry]
        H += out.ph_frac[entry] * out.SS_vec[entry].enthalpy # Xmol_phase * phase enthalpy
    end
    # Pure Phases
    if length(out.PP_vec) > 0
        for entry in 1:length(out.PP_vec)
            XH2Osolid += out.ph_frac[entry + length(out.SS_vec)] * out.PP_vec[entry].Comp[6] # Xmol_phase * H₂O content
            solid_n += out.ph_frac[entry + length(out.SS_vec)]
            H += out.ph_frac[entry + length(out.SS_vec)] * out.PP_vec[entry].enthalpy # Xmol_phase * phase enthalpy
        end
    end
    H /= solid_n
    XH2Osolid /= solid_n
    return XH2Osolid, H
end

# ====================================
# ========== H₂O EoS solver ==========
# ====================================
    # Pitzer and Sterner (1994) coefficients and EoS solver
    cᵢⱼ = [
        0.0               0.0             0.24657688e+6        0.51359951e+2        0.0                 0.0
        0.0               0.0             0.58638965e+0       -0.28646939e-2        0.31375577e-4       0.0
        0.0               0.0            -0.62783840e+1        0.14791599e-1        0.35779579e-3       0.15432925e-7
        0.0               0.0             0.0                 -0.42719875e+0       -0.16325155e-4       0.0
        0.0               0.0             0.56654978e+4       -0.16580167e+2        0.76560762e-1       0.0
        0.0               0.0             0.0                  0.10917883e+0        0.0                 0.0
        0.38878656e+13   -0.13494878e+9   0.30916564e+6        0.75591105e+1        0.0                 0.0
        0.0               0.0            -0.65537898e+5        0.18810675e+3        0.0                 0.0
        -0.14182435e+14   0.18165390e+9  -0.19769068e+6       -0.23530318e+2        0.0                 0.0
        0.0               0.0             0.92093375e+5        0.12246777e+3        0.0                 0.0
    ]
    get_cᵢ(cᵢⱼ, Tᵢ) = cᵢⱼ * [Tᵢ^(-4) Tᵢ^(-2) Tᵢ^(-1) 1.0 Tᵢ Tᵢ^2]'
    Pitzer_and_Sterner_EoS(cᵢ, ρ, Tᵢ) = (ρ + cᵢ[1]*ρ^2 - ρ^2*((cᵢ[3] + 2cᵢ[4]*ρ + 3cᵢ[5]*ρ^2 + 4cᵢ[6]*ρ^3) / (cᵢ[2] + cᵢ[3]*ρ + cᵢ[4]*ρ^2 + cᵢ[5]*ρ^3 + cᵢ[6]*ρ^4)^2) + cᵢ[7]*ρ^2*exp(-cᵢ[8]*ρ) + cᵢ[9]*ρ^2*exp(-cᵢ[10]*ρ))*(R*Tᵢ)

    # Iterative solver
    function PVT(P_target, Tᵢ, cᵢⱼ; verbose=false, volume=false)
        # Initial guess
        Ubound = 1e3; Lbound = 1e-3
        ρ = 0.5(Ubound+Lbound) # mol/cm³
        # Get coefficients
        cᵢ = get_cᵢ(cᵢⱼ, Tᵢ)
        # Iterative solve
        P_guess = 0.0
        error = 1e5
        for n in 1:70
            # Pitzer and Sterner EoS
            P_guess = Pitzer_and_Sterner_EoS(cᵢ, ρ, Tᵢ) # Pa
            # Update guess
            (P_guess > P_target) ? (Ubound = ρ) : (Lbound = ρ)
            ρ = 0.5(Lbound+Ubound)
            # Convergence check
            error = abs(P_guess - P_target)
        end
        (verbose) && println(" Pitzer and Sterner (1994) at $Tᵢ K and $(round(1e-5P_target, digits=4)) bar ($(round(1e-9P_target, digits=4)) GPa) \t:\t ρ = $(MM[end]*ρ) g/cm³ \t|\t V = $(1/ρ) cm³/mol \t (P error = $error Pa)")
        volume ? (return 1/ρ/MM[end]) : (return MM[end]*ρ)# in cm³/g or g/cm³
    end

    # Holland and Powell (2018) H₂O caloric parameterization
    const HP11 = (-241810.0, -298100.7, 188.8, 40.1, 0.8656e-2, 487500, -251.2) # Database Parameters (H⁰, G⁰, S⁰, c1, c2, c3, c5)
    # Reference state (1 bar, 25 °C)
    const Pᵣ = 1e5
    const Tᵣ = 298.15

    # Auxilliaries
        ∫CₚdT(c1, c3, c5, T, Tᵣ) = c1*(T-Tᵣ) - c3*(1/T - 1/Tᵣ) + 2c5*(√T - √Tᵣ) # J mol⁻¹
        ∫Cₚ_TdT(c1, c3, c5, T, Tᵣ) = c1*log(T/Tᵣ) - 0.5c3*(1/(T^2) - 1/(Tᵣ^2)) + 2c5*(1/√T - 1/√Tᵣ) # J mol⁻¹ K⁻¹
        ∫V_dP(P, Pᵣ, T, cᵢⱼ) = 1e-6(quadgk(x -> PVT(x, T, cᵢⱼ; volume=true), Pᵣ, P)[1]) # cm³ mol⁻¹ Pa -> (1e-6m³) mol⁻¹ Pa -> J mol⁻¹
        #∫V_dP(P, Pᵣ, T, cᵢⱼ; ns=1000) = @trapz LinRange(Pᵣ, P, ns) x ( 1e-6PVT(x, T, cᵢⱼ; volume=true) ) # cm³ mol⁻¹ Pa -> (1e-6m³) mol⁻¹ Pa -> J mol⁻¹

    # H₂O molar enthalpy, entropy, and gibbs energy
        H⁰ₚᵣ(P, T, Tᵣ) = HP11[1] + ∫CₚdT(HP11[4], HP11[6], HP11[7], T, Tᵣ) # H⁰(P,T) = H⁰(Pᵣ, T) + ∫CₚdT   (J mol⁻¹)
        S⁰ₚᵣ(P,T,Tᵣ) = HP11[3] + ∫Cₚ_TdT(HP11[4], HP11[6], HP11[7], T, Tᵣ) # S⁰(P,T) = S⁰(Pᵣ, T) + ∫CₚT⁻¹dT   (J mol⁻¹ K⁻¹)
    # H₂O molar gibbs solver
    # Note: P enters as GPa.
        g⁰(P,Pᵣ,T,cᵢⱼ) = (H⁰ₚᵣ(P, T, Tᵣ) - T*S⁰ₚᵣ(P, T, Tᵣ)) + ∫V_dP(1e9P, Pᵣ, T, cᵢⱼ) # J mol⁻¹

    # Numerical derivative of G⁰/T with respect to T at constant P
        δg⁰T⁻¹_δT(Pᵢ,Pᵣ,Tᵢ,cᵢⱼ) = ForwardDiff.derivative(Tᵢ -> g⁰(Pᵢ,Pᵣ,Tᵢ,cᵢⱼ)/Tᵢ, Tᵢ)

# Venusian Compositions
C = Dict(
    # Oxide list
    "oxides" => ["SiO2", "MgO", "FeO", "Fe2O3", "CaO", "Al2O3", "TiO2", "MnO", "K2O", "Na2O", "H2O"],
    # Compositions in wt%
    "VEGA"          => [50.8,  12.8, 7.7, 0.0, 8.3, 17.8, 0.22, 0.16, 0.11, 2.0, 0.0],          # Kargel et al. (1993), VEGA
    "VEGA-2"        => [45.6, 11.5, 7.7, 0.0, 7.5, 16, 0.2, 0.14, 0.1, 2.0, 0.0],               # Surkov et al. (1986), VEGA-2
    "VBasalt"       => [48.7, 8.1, 7.92, 0.98, 10.3, 17.9, 1.3, 8.1, 0.2, 2.4, 0.0],            # Semprich et al. (2024), Venusian Basalt
    "Venera-13"     => [45.1, 11.4, 9.3, 0.0, 7.1, 15.8, 1.6, 0.2, 4.0, 2.0, 0.0],              # Filiberto et al. (2014); Fegley (2003), Venera-13
    "Venera-14"     => [48.7, 8.1, 8.8, 0.0, 10.3, 17.9, 1.3, 0.2, 0.2, 2.4, 0.0],              # Filiberto et al. (2014); Fegley (2003), Venera-14
    "BSE"           => [45.0, 37.8, 8.05, 0.0, 3.55, 4.45, 0.201, 0.135, 0.029, 0.36, 0.0],     # McDonough and Sun (1995), Bulk Silicate Earth
)

# ============================
# ========== Inputs ==========
# ============================
# Add H₂O (mol%) → if vector, general solution is obtained by solving stacked P-T-H grids
H₂O_wt  = LinRange(0.0, 5.0, 10)
# H₂O_wt = 3.0
# Bulk composition choice
oX       = C["BSE"]
# Other parameters
Prange, Trange  = (0.0080, 1.0), (740, 1300)  # P (GPa) and T (K) ranges
nP, nT  = 80, 50  # Number of P-T points
generate_plots   = true
    plot_freeH2O = true
# ============================
# ====== Filled Outputs ======
# ============================
H₂O_mol = Hwt2mol(H₂O_wt, oX); nH = length(H₂O_mol)
nPnT, nPnTnH  = nP*nT, nP*nT*nH
output = Dict(  "Sys-ρ(kg/m³)"      => zeros(nPnTnH), "Sys-H⁰(kJ/mol)"     => zeros(nPnTnH), 
                "Solid-ρ(kg/m³)"    => zeros(nPnTnH), "Solid-H⁰(kJ/mol)"   => zeros(nPnTnH),
                "Liquid-ρ(kg/m³)"   => zeros(nPnTnH), "Liquid-H⁰(kJ/mol)"  => zeros(nPnTnH),
                "XH2Osolid"         => zeros(nPnTnH), "Vtotal(cm³/mol)"    => zeros(nPnTnH),
                "α_solid(1/K)"      => zeros(nPnTnH), "β_solid(1/GPa)"     => zeros(nPnTnH),
                "α_liquid(1/K)"     => zeros(nPnTnH), "β_liquid(1/GPa)"    => zeros(nPnTnH),
                "cptot(kJ/K/mol)"    => zeros(nPnTnH))
flvec = zeros(nPnTnH);
phi = zeros(nPnTnH);
# ============================
# ========= Algorithm ========
# ============================
function solve!(nP, nT, nH, nPnT, nPnTnH, Prange, Trange, H₂O_wt, oX, C, output, flvec, phi, generate_plots)
    attach_oX = zeros(nPnTnH); attach_oX[1:length(oX),end] .= wt2mol(oX)
    X    = [copy(x) for _ in H₂O_wt for x in map(Vector, eachrow(repeat(oX', outer=nPnT)))] # Vectorize compositions ensuring different pointers
    [X[i][end] = H₂O_wt[div(i-1, nPnT)+1] for i in eachindex(X)] # Assign H₂O contents
    X .= wt2mol.(X)
    p, t = LinRange(Prange[1], Prange[2], nP), LinRange(Trange[1], Trange[2], nT)  # Pressure in GPa --> Venusian surface sits at ~92 bar.
    if nH > 1
        P, T, H = repeat(p, nT*nH), repeat(repeat(t, inner=nP), outer=nH), repeat(H₂O_mol, inner=nP*nT) # Full P-T(-H₂O) grid
    else
        P, T = repeat(p, nT), repeat(t, inner=nP) # Full P-T grid
    end
    M = nH == 1 ? ["T(K)" "P(GPa)"] : ["T(K)" "P(GPa)" "H₂O_mol"]
    # Checks
    @assert minimum(T) >= 350 "T range must be >= 350 K"
    @assert minimum(P) >= 1e-4 "P range must be >= 1e-4 GPa (1 bar)"

    # Steps
    dP, dT, dH = step(p), step(t), nH>1 ? (H₂O_mol[2]-H₂O_mol[1]) : 0.0

    # # Parallel Minimization
    data    = Initialize_MAGEMin("ume", verbose=false);
    out = multi_point_minimization(10P, T.-273.15, data, X=X, Xoxides=C["oxides"], sys_in="mol", name_solvus=true) # takes kbar and °C
    Finalize_MAGEMin(data)

    # # Extract properties and assess pure H₂O
    Threads.@threads for i in eachindex(P)
            # out[].rho -> density in kg/m³
            # hₚₕ = out[].enthalpy -> single-element vector containing the enthalpy (per kJ/kg or kJ/mol, depending on 'sys_in') of the system

            # XH2Osolid and solid enthalpy
            XH2Os, sH = XH2Osolid_sH(out[i])
            # Liquid phase
            output["Liquid-ρ(kg/m³)"][i] = 1e3PVT(1e9*P[i], T[i], cᵢⱼ; verbose=false) # (kg m⁻³)
            output["Liquid-H⁰(kJ/mol)"][i] = -1e-3(T[i]^2)δg⁰T⁻¹_δT(P[i] ,Pᵣ,T[i],cᵢⱼ) # (kJ/mol)
            # Solid phase
            output["Solid-ρ(kg/m³)"][i] = out[i].rho_S #  (kg m⁻³)
            output["Solid-H⁰(kJ/mol)"][i] = sH # (kJ/mol)
            # Minimized data (System)
            output["Sys-ρ(kg/m³)"][i] = output["Solid-ρ(kg/m³)"][i]*out[i].frac_S + output["Liquid-ρ(kg/m³)"][i]*out[i].frac_F #  (kg m⁻³)
            output["Sys-H⁰(kJ/mol)"][i] = output["Solid-H⁰(kJ/mol)"][i]*out[i].frac_S + output["Liquid-H⁰(kJ/mol)"][i]*out[i].frac_F # (kJ/mol)
            # Minimized volume
            output["Vtotal(cm³/mol)"][i] = out[i].V # cm³ mol⁻¹
            # Total cₚ
            output["cptot(kJ/K/mol)"][i] = out[i].s_cp[1] # kJ K⁻¹ kg⁻¹
            # Minimized XH2Osolid
            output["XH2Osolid"][i] = XH2Os # mol fraction
            # Fluid fraction
            flvec[i] = out[i].frac_F
            phi[i] = (output["Sys-ρ(kg/m³)"][i] - output["Solid-ρ(kg/m³)"][i])/(output["Liquid-ρ(kg/m³)"][i] - output["Solid-ρ(kg/m³)"][i])
    end

    # Dummy reshapes
    ρs = reshape(output["Solid-ρ(kg/m³)"], nP, nT, nH)
    ρl = reshape(output["Liquid-ρ(kg/m³)"], nP, nT, nH)
    αs, αl, βs, βl = zeros(nP,nT,nH), zeros(nP,nT,nH), zeros(nP,nT,nH), zeros(nP,nT,nH)

    # Compute α and β
    Threads.@threads for iH in eachindex(H₂O_mol)
        for iT in eachindex(t)
            for iP in eachindex(p)
                # αs
                (iT==1)   &&  (αs[iP,iT,iH] = -(1/ρs[iP,iT,iH]) * (ρs[iP,iT+1,iH]-ρs[iP,iT,iH])/dT) # left bc
                (iT==nT)  &&  (αs[iP,iT,iH] = -(1/ρs[iP,iT,iH]) * (ρs[iP,iT,iH]-ρs[iP,iT-1,iH])/dT) # right bc
                (iT>1 && iT<nT)  &&  (αs[iP,iT,iH] = -(1/ρs[iP,iT,iH]) * (ρs[iP,iT+1,iH]-ρs[iP,iT-1,iH])/(2dT))
                # αl
                (iT==1)   &&  (αl[iP,iT,iH] = -(1/ρl[iP,iT,iH]) * (ρl[iP,iT+1,iH]-ρl[iP,iT,iH])/dT) # left bc
                (iT==nT)  &&  (αl[iP,iT,iH] = -(1/ρl[iP,iT,iH]) * (ρl[iP,iT,iH]-ρl[iP,iT-1,iH])/dT) # right bc
                (iT>1 && iT<nT)  &&  (αl[iP,iT,iH] = -(1/ρl[iP,iT,iH]) * (ρl[iP,iT+1,iH]-ρl[iP,iT-1,iH])/(2dT))
                # βs
                (iP==1)   &&  (βs[iP,iT,iH] = (1/ρs[iP,iT,iH]) * (ρs[iP+1,iT,iH]-ρs[iP,iT,iH])/dP) # left bc
                (iP==nP)  &&  (βs[iP,iT,iH] = (1/ρs[iP,iT,iH]) * (ρs[iP,iT,iH]-ρs[iP-1,iT,iH])/dP) # right bc
                (iP>1 && iP<nP)  &&  (βs[iP,iT,iH] = (1/ρs[iP,iT,iH]) * (ρs[iP+1,iT,iH]-ρs[iP-1,iT,iH])/(2dP))
                # βl
                (iP==1)   &&  (βl[iP,iT,iH] = (1/ρl[iP,iT,iH]) * (ρl[iP+1,iT,iH]-ρl[iP,iT,iH])/dP) # left bc
                (iP==nP)  &&  (βl[iP,iT,iH] = (1/ρl[iP,iT,iH]) * (ρl[iP,iT,iH]-ρl[iP-1,iT,iH])/dP) # right bc
                (iP>1 && iP<nP)  &&  (βl[iP,iT,iH] = (1/ρl[iP,iT,iH]) * (ρl[iP+1,iT,iH]-ρl[iP-1,iT,iH])/(2dP))
            end
        end
    end
    output["α_solid(1/K)"]    .= vec(αs)
    output["α_liquid(1/K)"]   .= vec(αl)
    output["β_solid(1/GPa)"]  .= vec(βs)
    output["β_liquid(1/GPa)"] .= vec(βl)

    # Field Rounding
    P .= round.(P, digits=4)
    T .= round.(T, digits=4)
    (nH>1) && (H .= round.(H, digits=4))
    output["Sys-ρ(kg/m³)"]       .= round.(output["Sys-ρ(kg/m³)"], digits=4)
    output["Sys-H⁰(kJ/mol)"]     .= round.(output["Sys-H⁰(kJ/mol)"], digits=4)
    output["Solid-ρ(kg/m³)"]     .= round.(output["Solid-ρ(kg/m³)"], digits=4)
    output["Solid-H⁰(kJ/mol)"]   .= round.(output["Solid-H⁰(kJ/mol)"], digits=4)
    output["Liquid-ρ(kg/m³)"]    .= round.(output["Liquid-ρ(kg/m³)"], digits=4)
    output["Liquid-H⁰(kJ/mol)"]  .= round.(output["Liquid-H⁰(kJ/mol)"], digits=4)
    output["XH2Osolid"]          .= round.(output["XH2Osolid"], digits=6)
    output["Vtotal(cm³/mol)"]    .= round.(output["Vtotal(cm³/mol)"], digits=6)
    output["α_solid(1/K)"]       .= round.(output["α_solid(1/K)"], digits=8)
    output["β_solid(1/GPa)"]     .= round.(output["β_solid(1/GPa)"], digits=8)
    output["α_liquid(1/K)"]      .= round.(output["α_liquid(1/K)"], digits=8)
    output["β_liquid(1/GPa)"]    .= round.(output["β_liquid(1/GPa)"], digits=8)

    # # Phase iterator + matrix constructor
    !isdir("./+op") && mkdir("./+op")
    if nH == 1
        # Check observed phases
        observed_phases = Dict{String, Array{Float64, 2}}()
        for n in eachindex(out)
            # Vectorized to mesh indexing
            j = Int(ceil(n/nP)); i = n - (j-1)*nP;
            # Phase iterator
            for (idx, phase) in enumerate(out[n].ph)
                # Initialise phase if not present
                !haskey(observed_phases, phase) && (observed_phases[phase] = zeros(Float64, nP, nT))
                # Assign solution phase value
                observed_phases[phase][i, j] = out[n].ph_frac_wt[idx]
            end
        end
        println("Observed phases in the P-T grid: ", collect(keys(observed_phases)))

        # Construct vectorized data matrix
        M = vcat(M, [T P]); 
        for (key, values) in output; M = hcat(M, vcat(key, values)); end
        M = hcat(M, vcat("X(mol_frac)", attach_oX))
    else
        # Check observed phases
        observed_phases = String[]
        for n in eachindex(out)
            for phase in out[n].ph
                !(phase in observed_phases) && push!(observed_phases, phase)
            end
        end
        println("Observed phases in the P-T-H grid : ", observed_phases)

        # Construct vectorized data matrix
        M = vcat(M, [T P 1e-2H]); 
        for (key, values) in output; M = hcat(M, vcat(key, values)); end
        M = hcat(M, vcat("X(mol_frac)", attach_oX))
    end

    # Write to file
    open("./+op/table.dat", "w") do io; writedlm(io, M, '\t'); end

    # # Plots
    if generate_plots

        # Plot 
        fsz = (1200, 800)
        interp = true
        cmap = [:OrRd, Reverse(:davos), :algae, :PuRd, :Purples, :PuOr, :diff, :PRGn, :BrBg, :vik100, :vanimo, :tofino, :lisbon, :glasgow, :cork, :batlow, :bamako]

        if nH == 1

            xlabsz, ylabsz, titlesz, xticklabsz, yticklabsz, xticksz, yticksz = 20, 20, 22, 16, 16, 12, 12
            cmap = :vik100

            # Mineral phase map
            fig = Figure(size=fsz); i=1; j=1; ncols=4
            for (ph, map) in observed_phases
                ax = Axis(fig[i, j], title=ph, ylabel=L"Pressure\;[GPa]", xlabel=L"Temperature\;[K]", xlabelsize=xlabsz, ylabelsize=ylabsz,
                            xticklabelsize=xticklabsz, yticklabelsize=yticklabsz, xticksize=xticksz, yticksize=yticksz, titlesize=titlesz, yreversed=true)
                hm = heatmap!(ax, t, p, map', colormap=Reverse(:davos50), interpolate=true)
                Colorbar(fig[i, j+1], hm; labelsize=20, ticklabelsize=15, vertical=true, labelpadding=13)
                j+=2; (j==2ncols+1) && (i+=1; j=1)
            end
            GLMakie.save("./Minmap.png", fig)

            # System density
            fig = Figure(size = fsz)
            ax1 = Axis(fig[1, 1], xlabel=L"T\;(K)", ylabel=L"P\;(GPa)", title="System ρ (kg/m³)", yreversed=true, xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz)
            hm = heatmap!(ax1, unique(T), unique(P), reshape(output["Sys-ρ(kg/m³)"], nP, nT)', colormap = cmap, interpolate=interp)
            Colorbar(fig[1, 2], hm)

            # System enthalpy
            ax2 = Axis(fig[2, 1], xlabel=L"T\;(K)", ylabel=L"P\;(GPa)", title="System Enthalpy (kJ/mol)", yreversed=true, xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz)
            hm = heatmap!(ax2, unique(T), unique(P), reshape(output["Sys-H⁰(kJ/mol)"], nP, nT)', colormap = cmap, interpolate=interp)
            Colorbar(fig[2, 2], hm)

            # H₂O density
            ax3 = Axis(fig[1, 3], xlabel=L"T\;(K)", ylabel=L"P\;(GPa)", title="H₂O/Liquid ρ (kg/m³)", yreversed=true, xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz)
            hm = heatmap!(ax3, unique(T), unique(P), reshape(output["Liquid-ρ(kg/m³)"], nP, nT)', colormap = cmap, interpolate=interp)
            Colorbar(fig[1, 4], hm)

            # H₂O enthalpy
            ax4 = Axis(fig[2, 3], xlabel=L"T\;(K)", ylabel=L"P\;(GPa)", title="H₂O/Liquid Enthalpy (kJ/mol)", yreversed=true, xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz)
            hm = heatmap!(ax4, unique(T), unique(P), reshape(output["Liquid-H⁰(kJ/mol)"], nP, nT)', colormap = cmap, interpolate=interp)
            Colorbar(fig[2, 4], hm)

            # Solid density
            ax3 = Axis(fig[1, 5], xlabel=L"T\;(K)", ylabel=L"P\;(GPa)", title="Solid ρ (kg/m³)", yreversed=true, xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz)
            hm = heatmap!(ax3, unique(T), unique(P), reshape(output["Solid-ρ(kg/m³)"], nP, nT)', colormap = cmap, interpolate=interp)
            Colorbar(fig[1, 6], hm)

            # Solid enthalpy
            ax4 = Axis(fig[2, 5], xlabel=L"T\;(K)", ylabel=L"P\;(GPa)", title="Solid Enthalpy (kJ/mol)", yreversed=true, xlabelsize=xlabsz, ylabelsize=ylabsz, titlesize=titlesz)
            hm = heatmap!(ax4, unique(T), unique(P), reshape(output["Solid-H⁰(kJ/mol)"], nP, nT)', colormap = cmap, interpolate=interp)
            Colorbar(fig[2, 6], hm)

            # Save
            GLMakie.save("./Hρ.png", fig)

            # display
            display(fig)

        else

            xlabsz, ylabsz, titlesz, xticklabsz, yticklabsz, xticksz, yticksz, zlabsz, zticklabsz, zticksz = 15, 15, 22, 10, 10, 12, 10, 15, 12, 10
            wf_subsample = 2

            if plot_freeH2O

                fig = Figure(size = fsz)
                ax = Axis3(fig[1, 1], xlabel="T", ylabel="P", zlabel="H", xlabelsize=xlabsz, ylabelsize=ylabsz, zlabelsize=zlabsz, titlesize=titlesz)
                scatter!(ax, M[2:end,1], M[2:end,2], M[2:end,3], markersize=10, alpha=1.0, colormap=:Blues, color=flvec)
                display(fig)

            else
                # System density 3D surface
                fig = Figure(size = fsz)
                ax = Axis3(fig[2, 1], xlabel=L"T\;(\mathrm{K})", ylabel=L"P\;(\mathrm{GPa})", zlabel=L"System\;\rho\;[\mathrm{kg/m^3}]", xlabelsize=xlabsz, ylabelsize=ylabsz, zlabelsize=zlabsz, titlesize=titlesz)
                smap = reshape(output["Sys-ρ(kg/m³)"], nP, nT, nH)
                for (i, h₂o) in enumerate(H₂O_wt)
                    surface!(ax, unique(T), unique(P), smap[:,:,i]', label = "$h₂o", colormap=cmap[i])
                    wireframe!(ax, unique(T)[1:wf_subsample:end], unique(P)[1:wf_subsample:end], smap[:,:,i][1:wf_subsample:end,1:wf_subsample:end]', color = :black, linewidth = 0.5, alpha = 0.3)
                end

                # Liquid density 3D surface
                ax = Axis3(fig[2, 2], xlabel=L"T\;(\mathrm{K})", ylabel=L"P\;(\mathrm{GPa})", zlabel=L"Liquid\;\rho\;[\mathrm{kg/m^3}]", xlabelsize=xlabsz, ylabelsize=ylabsz, zlabelsize=zlabsz, titlesize=titlesz)
                smap = reshape(output["Liquid-ρ(kg/m³)"], nP, nT, nH)
                for (i, h₂o) in enumerate(H₂O_wt)
                    surface!(ax, unique(T), unique(P), smap[:,:,i]', label = "$h₂o", colormap=cmap[i])
                    wireframe!(ax, unique(T)[1:wf_subsample:end], unique(P)[1:wf_subsample:end], smap[:,:,i][1:wf_subsample:end,1:wf_subsample:end]', color = :black, linewidth = 0.5, alpha = 0.3)
                end

                # Solid density 3D surface
                ax = Axis3(fig[2, 3], xlabel=L"T\;(\mathrm{K})", ylabel=L"P\;(\mathrm{GPa})", zlabel=L"Solid\;\rho\;[\mathrm{kg/m^3}]", xlabelsize=xlabsz, ylabelsize=ylabsz, zlabelsize=zlabsz, titlesize=titlesz)
                smap = reshape(output["Solid-ρ(kg/m³)"], nP, nT, nH)
                for (i, h₂o) in enumerate(H₂O_wt)
                    surface!(ax, unique(T), unique(P), smap[:,:,i]', label = "$h₂o", colormap=cmap[i])
                    wireframe!(ax, unique(T)[1:wf_subsample:end], unique(P)[1:wf_subsample:end], smap[:,:,i][1:wf_subsample:end,1:wf_subsample:end]', color = :black, linewidth = 0.5, alpha = 0.3)
                end

                # System Enthalpy 3D surface
                ax = Axis3(fig[3, 1], xlabel=L"T\;(\mathrm{K})", ylabel=L"P\;(\mathrm{GPa})", zlabel=L"System\;H⁰\;[\mathrm{kJ/mol}]", xlabelsize=xlabsz, ylabelsize=ylabsz, zlabelsize=zlabsz, titlesize=titlesz)
                smap = reshape(output["Sys-H⁰(kJ/mol)"], nP, nT, nH)
                for (i, h₂o) in enumerate(H₂O_wt)
                    surface!(ax, unique(T), unique(P), smap[:,:,i]', label = "$h₂o", colormap=cmap[i])
                    wireframe!(ax, unique(T)[1:wf_subsample:end], unique(P)[1:wf_subsample:end], smap[:,:,i][1:wf_subsample:end,1:wf_subsample:end]', color = :black, linewidth = 0.5, alpha = 0.3)
                end

                # Liquid density 3D surface
                ax = Axis3(fig[3, 2], xlabel=L"T\;(\mathrm{K})", ylabel=L"P\;(\mathrm{GPa})", zlabel=L"Liquid\;H⁰\;[\mathrm{kJ/mol}]", xlabelsize=xlabsz, ylabelsize=ylabsz, zlabelsize=zlabsz, titlesize=titlesz)
                smap = reshape(output["Liquid-H⁰(kJ/mol)"], nP, nT, nH)
                for (i, h₂o) in enumerate(H₂O_wt)
                    surface!(ax, unique(T), unique(P), smap[:,:,i]', label = "$h₂o", colormap=cmap[i])
                    wireframe!(ax, unique(T)[1:wf_subsample:end], unique(P)[1:wf_subsample:end], smap[:,:,i][1:wf_subsample:end,1:wf_subsample:end]', color = :black, linewidth = 0.5, alpha = 0.3)
                end

                # Solid density 3D surface
                ax = Axis3(fig[3, 3], xlabel=L"T\;(\mathrm{K})", ylabel=L"P\;(\mathrm{GPa})", zlabel=L"Solid\;H⁰\;[\mathrm{kJ/mol}]", xlabelsize=xlabsz, ylabelsize=ylabsz, zlabelsize=zlabsz, titlesize=titlesz)
                smap = reshape(output["Solid-H⁰(kJ/mol)"], nP, nT, nH)
                for (i, h₂o) in enumerate(H₂O_wt)
                    surface!(ax, unique(T), unique(P), smap[:,:,i]', label = "$h₂o", colormap=cmap[i])
                    wireframe!(ax, unique(T)[1:wf_subsample:end], unique(P)[1:wf_subsample:end], smap[:,:,i][1:wf_subsample:end,1:wf_subsample:end]', color = :black, linewidth = 0.5, alpha = 0.3)
                end

                # Legend
                Legend(fig[1, 1:3], ax, "--- H₂O [wt%] ---"; labelsize=20, ticklabelsize=12, orientation=:horizontal, framevisible=true, titlesize=20)
            end

            # display
            display(fig)
        end
    end

    return out
end
out = solve!(nP, nT, nH, nPnT, nPnTnH, Prange, Trange, H₂O_wt, oX, C, output, flvec, phi, generate_plots)