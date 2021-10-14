function T(λ,Q0,tspan)
    D = Derivative()
    B = Evaluation(tspan[1])
    A = [B               0;
         B*D             0;
         0               B;
         0               B*D;
         -D^2-(Q0+λ)*I   0;
         -I              -D^2-(Q0+λ)*I];
    u,udλ = A\[1;0;0;0;0;0]
    v,vdλ = A\[0;1;0;0;0;0]
    return [u(tspan[2]) v(tspan[2]);
    u'(tspan[2]) v'(tspan[2]);
    udλ(tspan[2]) vdλ(tspan[2]);
    udλ'(tspan[2]) vdλ'(tspan[2])]
end

function Δ(λ,Q0,tspan)
    Ts = T(λ,Q0,tspan)
   (Ts[1,1] + Ts[2,2])/2
end

function Δ2(λ,Q0,tspan)
    Ts = T(λ,Q0,tspan)
   ((Ts[1,1] + Ts[2,2])/2)^2
end

function Δ_newton_step(λ,σ,Q0,tspan)
    Ts = T(λ,Q0,tspan)
    δ = (Ts[1,1] + Ts[2,2] - 2*σ)/(Ts[3,1] + Ts[4,2])
    (λ - δ,δ)
end

function Δ2_newton_step(λ,Q0,tspan)
    Ts = T(λ,Q0,tspan)
    Δ = (Ts[1,1] + Ts[2,2])/2
    dΔ = (Ts[3,1] + Ts[4,2])/2
    δ = (Δ^2 - 1)/(2*Δ*dΔ)
    (λ - δ,δ)
end

function Δroot(λ,σ,Q0,tspan)
    λ0 = λ
    for i = 1:100
        (λ0,δ) = Δ_newton_step(λ0,σ,Q0,tspan)
        if abs(δ) < 1e-15
            break
        end
        if i == 100
            @warn "Max iter"
        end
    end
    λ0
end

function Δ2root(λ,Q0,tspan)
    λ0 = λ
    δ = 0.0
    for i = 1:200
        (λ0,δ) = Δ2_newton_step(λ0,Q0,tspan)
        if abs(δ) < 1e-14
            break
        end
        if i == 100
            @warn "Max iter"
        end
    end
    (λ0,δ)
end

function T12(λ,Q0,tspan)
    Ts = T(λ,Q0,tspan)
    Ts[1,2]
end

function T12_newton_step(λ,Q0,tspan)
    Ts = T(λ,Q0,tspan)
    δ = Ts[1,2]/Ts[3,2]
    (λ - δ,δ)
end

function T12root(λ,Q0,tspan)
    λ0 = λ
    δ = 0.0
    for i = 1:100
        (λ0,δ) = T12_newton_step(λ0,Q0,tspan)
        if abs(δ) < 1e-14
            break
        end
        if i == 100
            @warn "Max iter"
        end
    end
    (λ0,δ)
end

sqz = z -> sqrt(z-1 |> complex)*sqrt(1+z |> complex)

function find_sheet(λ,Q0,tspan)
    TT = T(λ+1im*1e-14,Q0,tspan)
    dd = (TT[1,1] + TT[2,2])/2
    sqdd = sqz(dd)
    #return -sign(imag(dd^2))
    #println(imag(dd^2) |> sign)
    s1 = sqdd + 1/2*(TT[2,2]-TT[1,1]) |> abs
    s2 = -sqdd + 1/2*(TT[2,2]-TT[1,1]) |> abs
    #println((s1,s2))
    if s1 > s2 # pole on sheet 1
        σ = 1.0
    else
        σ = -1.0 # pole on sheet 2
    end
    return σ
end

function DirichletSpectrum(q0,L,n,tol)
    xspan = (0.0,L)
    Q0 = Fun(q0,xspan[1]..xspan[2])
    sp = Q0.space
    BC = Dirichlet(sp)
    D = Derivative()
    L = -D^2 - Q0
    λDir, v = ApproxFun.eigs(BC,L,n, tolerance = tol);
    λDir |> sort
end

function PeriodicSpectrum(q0,L,n,tol)
    d = ApproxFun.PeriodicSegment(0.0,2*L)
    Q0per = Fun(q0,d)
    D = Derivative(Q0per.space)
    L = -D^2 - Q0per
    λPer, v = ApproxFun.eigs(L,n, tolerance = tol);
    λPer |> sort
end

function ScatteringData(q0,L,n,tol,k::Integer)
    xspan = (0.0,L)
    Q0 = Fun(q0,xspan[1]..xspan[2])

    if 2k>n
        @warn "Need to increase n"
    end
    λPer = PeriodicSpectrum(q0,L,n,tol)
    λDir = DirichletSpectrum(q0,L,n,tol)
    σs = map(y -> find_sheet(y,Q0,xspan), λDir[1:k]);
    α1 = λPer[1]
    α = λPer[3:2:2k+1] .- α1
    β = λPer[2:2:2k] .- α1
    gaps = hcat(β,α)
    zs = hcat(λDir[1:k] .- α1, σs)
    (gaps,zs,α1)
end

function ScatteringData(q0,L,n,tol,trunctol::Float64)
    xspan = (0.0,L)
    Q0 = Fun(q0,xspan[1]..xspan[2])
    λPer = PeriodicSpectrum(q0,L,n,tol)
    α1 = λPer[1]
    α = λPer[3:2:end] .- α1
    β = λPer[2:2:end] .- α1
    m = min(length(α),length(β))
    gaplen = - β[1:m] + α[1:m]
    k = 0
    for j = 1:m
        k = j-1
        if gaplen[j] < trunctol # TODO: Is this right?
            break
        end
    end
    #display(gaplen)

    λDir = DirichletSpectrum(q0,L,2*n,tol)
    σs = map(y -> find_sheet(y,Q0,xspan), λDir[1:k]);
    α1 = λPer[1]
    α = λPer[3:2:2k+1] .- α1
    β = λPer[2:2:2k] .- α1
    gaps = hcat(β,α)
    zs = hcat(λDir[1:k] .- α1, σs)
    (gaps,zs,α1)
end



function find_s(cΩx::Vector,cΩ0::Vector,L,gaps::Array,tol)
    Ωx = imag(cΩx); Ω0 = imag(cΩ0)
    metric = sqrt.(gaps[:,2] - gaps[:,1])
    dt = 0.001*pi/L
    t = 0.0
    T = 0.0
    v =  Ω0
    Δ = 1e-8
    for j = 1:1000
        v += Ωx*dt
        v = mod2pi.(v)
        t += dt
        d1 = minimum(v./metric)
        d2 = minimum(abs.(v .- pi)./metric)
        d3 = minimum(abs.(v .- 2pi)./metric)
        if min(d1,d2,d3) > Δ
            Δ = min(d1,d2,d3)
            T = t
            if Δ > tol
                return mod(T,L)
            end
        end
    end
    @warn "Failed to find shift"
    return mod(T,L)
end

function RefineHyperellipticSurface!(q0::Function,L,S::HyperellipticSurface,invtol,n,tol,k,m,new=true)
    s = find_s(S.Ωx,S.Ω0,L,S.gaps,invtol)
    gaps, zs, α1 = PeriodicKdV.ScatteringData(x -> q0(x-s),L,n,tol,k)
    Snew = HyperellipticSurface(gaps,zs,α1,m,new)
    S.Ω0 = Snew.Ωx*s+Snew.Ω0
end
