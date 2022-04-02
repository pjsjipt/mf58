
using Gridap
using GridapGmsh
using Gridap.Geometry


mutable struct ModelInfo
    "Reference velocity"
    Uref::Float64
    "Velocity on different surfaces"
    U::Vector{Float64}
    "Different material regions"
    regions::Vector{String}
    "Surface regions where convection occurs"
    hsurf::Vector{String}
    "Surface regions where heat flux is prescribed"
    qsurf::Vector{String}
    "Surface regions where temperature is prescribed"
    tsurf::Vector{String}
    "Regions to which the prescribed temperature belongs"
    surfreg::Dict{String,String}
    "Heat conductivity on different materials"
    k::Dict{String,Float64}
    "Convection coefficient on different surfaces"
    h::Dict{String,Float64}
    "Prescribed temperature on different surfaces"
    T::Dict{String,Float64}
    "Prescribed heat flux on different surfaces"
    q::Dict{String,Float64}
    "Heat conductivity of fluif"
    kₐ::Float64
    "Prandtl number of fluid"
    Prₐ::Float64
    "Kinematic viscosity of the fluid"
    νₐ::Float64
    "Fluid temperature"
    Tref::Float64
end

"""
`ModelInfo(u=10.0; k=Dict("glass"=>1.0, "steel"=>50.0), hsurf=["hglass", "hsteel"],
                   T=Dict("temp_glass"=>200.0, "temp_steel"=>200.0),
                   q=Dict{String,Float64}(), Lmm = [2.0, 0.5],
                   qsurf=Dict{String,Float64}(),
                   Prₐ=0.7, kₐ=25e-3, νₐ=1.6e-5, Tref=20.0)`

Specifies material parameters and boundary conditions for FEM analysis. 
The heat conduction model assumes that several regions with differing materials
can exist. There are 3 possible boundary conditions

 1. Prescribed temperature, Dirichlet BC
 2. Prescribed heat flux, Neumann BC
 3. Prescribed convection coefficient, mixed BC

Boundaries where no BC is prescribed automatically uses a prescribed 
0 heat flux. This BC can be used for symmetry.

In the case of convection, different surfaces can have different velocities.

For now, to calculate the convection coefficient, a correlation for circular 
cylinders is used:

      ``Nu = 0.51⋅Re^(½)⋅Pr^(0.37)``


"""
function ModelInfo(u=10.0; k=Dict("glass"=>1.0, "steel"=>50.0), hsurf=["hglass", "hsteel"],
                   T=Dict("temp_glass"=>200.0, "temp_steel"=>200.0),
                   surfreg=Dict("temp_glass"=>"glass", "temp_steel"=>"steel",
                                "hglass"=>"glass", "hsteel"=>"steel"),
                   q=Dict{String,Float64}(), Lmm = [2.0, 0.5],
                   qsurf=Dict{String,Float64}(),
                   Prₐ=0.707, kₐ=25.3e-3, νₐ=1.6e-5, Tref=20.0)
    nh = length(hsurf)
    nu = length(u)

    Uref = u[1]

    if nu==1
        U = fill(u[1], nh)
    elseif nu != nh
        error("Incompatible number of convection surfaces with velocities!")
    else
        U = [uu for uu in u]
    end
    

    if length(Lmm) != nh
        error("The number of dimensions should be the same as the number of convection regions!")
    end
    
    L = Lmm .* 1e-3
    Re = (U ./ νₐ) .* L
    Nu = (0.51*Prₐ^0.37) .* Re .^ 0.35 * 1.35 # .^ 0.5
    h1 = kₐ .* Nu ./ L

    h = Dict{String,Float64}()

    for i in 1:nh
        h[ hsurf[i] ] = h1[i]
    end
    return  ModelInfo(Uref, U,  collect(keys(k)), hsurf, collect(keys(q)), collect(keys(T)),
                      surfreg, k, h, T, q, kₐ, Prₐ, νₐ, Tref)

end

r(x) = x[2]

function fem_model(th::ModelInfo, model; order=2, dimension=2)

    degree = 2*order
    reffe = ReferenceFE(lagrangian, Float64, order)

    # Vamos montar os espaçoes de funções
    

    tsurf = th.tsurf
    uw = [th.T[s]-th.Tref for s in tsurf]
    V = TestFESpace(model, reffe, dirichlet_tags=tsurf)
    U = TrialFESpace(V, uw)

    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    hsurf = th.hsurf  # Boundaries with convection
    Γh = [BoundaryTriangulation(model, tags=hs) for hs in hsurf]
    dΓh = [Measure(Γ,degree) for Γ in Γh]

    qsurf = th.qsurf
    #Γq = [BoundaryTriangulation(model, tags=hs) for qs in qsurf]
    #dΓq = [Measure(Γ,degree) for Γ in Γq]
    
    Γt = [BoundaryTriangulation(model, tags=ts) for ts in tsurf]
    dΓt = [Measure(Γ,degree) for Γ in Γt]

    labels = get_face_labeling(model)
    tags = get_face_tag(labels, dimension)
    kdict = Dict{Int,Float64}()
    for r in th.regions
        tagidx = get_tag_from_name(labels, r)
        kdict[tagidx] = th.k[r]
    end
    
    function heatflux(∇u, tag)
        k = kdict[tag]
        return k * ∇u
    end
    f(x) = 0.0
    nh = length(hsurf)
    hh = zeros(nh)
    for i in 1:nh
        hh[i] = th.h[ hsurf[i] ]
    end
    
    function a(u,v)
        z = ∫( r*∇(v)⋅(heatflux∘(∇(u), tags) ) )*dΩ
        for i in 1:nh
            z = z + ∫(r*v*u*hh[i])*dΓh[i]
        end
        return z
    end
    #return 
    b(v) = ∫(r*v*f)*dΩ

    op = AffineFEOperator(a,b,U,V)

    ls = LUSolver()
    solver = LinearFESolver(ls)

    uh = solve(solver, op)

    return ( u=uh, Ω=(Ω,dΩ), Γh=(Γh,dΓh), Γt=(Γt,dΓt) )
    
    
end




function process_thermistor(th, fem)

    (;u, Ω, Γh, Γt) = fem
    # Calculate surface areas (m²)
    Ah = 2π .* [sum(∫(r)*dΓ) for dΓ in Γh[2]]
    At = 2π .* [sum(∫(r)*dΓ) for dΓ in Γt[2]]

    # Calculate the mean temperature on the surfaces
    Th = 2π .* [sum(∫( r*u ) *dΓ ) for dΓ in Γh[2]] ./ Ah
        
    # Heat flux from surfaces:
    # Convection surfaces:
    Qh = Float64[]
    nh = length(th.hsurf)
    for i = 1:nh
        hs = th.hsurf[i]
        h = th.h[hs]
        dΓ = Γh[2][i]
        push!(Qh, 4π*h * sum(∫( r*u ) * dΓ))
    end
    # Heat flux from prescribed temperature regions
    Qt = Float64[]
    nt = length(th.tsurf)
    for i in 1:nt
        s = th.tsurf[i]
        k = th.k[th.surfreg[s]]
        Γ  = Γt[1][i]
        dΓ = Γt[2][i]
        n⃗ = get_normal_vector(Γ)
        push!(Qt, -4π*k * sum( ∫( r*(∇(u)⋅n⃗ ) ) * dΓ ))
    end

    # Heat flux from convection surfaces calculated from heat conduction
    Qh1 = Float64[]
    for i in 1:nh
        s = th.tsurf[i]
        k = th.k[th.surfreg[s]]
        Γ  = Γh[1][i]
        dΓ = Γh[2][i]
        n⃗ = get_normal_vector(Γ)
        push!(Qh1, -4π*k * sum( ∫( r*(∇(u)⋅n⃗ ) ) * dΓ ))
    end
    
    return Th, Qh, Qt, Qh1, Ah, At

    
end

