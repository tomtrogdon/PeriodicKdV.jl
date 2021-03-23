include("./PeriodicKdV.jl")
using .PeriodicKdV
using LinearAlgebra, Plots, Printf
using ApproxFun

q0 = x -> cos(2.0*x)
L = pi*1.0

S = PeriodicKdV.HyperellipticSurface(q0,L,100,10000.0,100,1e-14,.4)




BA = BakerAkhiezerFunction(S,20,1e-14)

S.Ωx



out = PeriodicKdV.PeriodicSpectrum(q0,L,100,1.0e-14)













out = PeriodicKdV.ScatteringData(q0,L,100,1.0e-14,10)

out[1][:,1] - out[1][:,2]















function mScatteringData(q0,L,n,tol,trunctol::Float64)
    xspan = (0.0,L)
    Q0 = Fun(q0,xspan[1]..xspan[2])
    λPer = PeriodicKdV.PeriodicSpectrum(q0,L,n,tol)
    α1 = λPer[1]
    α = λPer[3:2:end] .- α1
    β = λPer[2:2:end] .- α1
    m = min(length(α),length(β))
    gaplen = - β[1:m] + α[1:m]
    k = 0
    for j = 1:m
        k = j
        if gaplen[j] < trunctol # TODO: Is this right?
            break
        end
    end
    display(k)

    λDir = PeriodicKdV.DirichletSpectrum(q0,L,2*n,tol)
    σs = map(y -> PeriodicKdV.find_sheet(y,Q0,xspan), λDir[1:k]);
    α1 = λPer[1]
    α = λPer[3:2:2k+1] .- α1
    β = λPer[2:2:2k] .- α1
    gaps = hcat(β,α)
    zs = hcat(λDir[1:k] .- α1, σs)
    (gaps,zs,α1)
end



out = mScatteringData(q0,L,100,1.0e-14,1.0e-13)




out = PeriodicKdV.DirichletSpectrum(q0,L,100,1.0e-14)

u(x) = u(x + L)
u(x) = -u(x + L)


out[1]








PeriodicKdV.ScatteringData(q0,L,100,1e-10,1e-10)





q0(0.3) - KdV(BA,0.3,0.0)
