function circle_points(a,r,n)
    z = .5:1.0:n-.5 |> Array
    z *= (2im*pi)/n
    z = exp.(z)
    r*z.+a
end

struct CircleFun
    a::ComplexF64
    r::Float64
    v::Vector{ComplexF64}
end

function CircleFun(f,a,r,n)
    CircleFun(a,r,f(circle_points(a,r,n)))
end

function Integrate(F::CircleFun)
    sum(F.v.*(circle_points(F.a,F.r,length(F.v)) .- F.a))*(2im*pi)/length(F.v)
end

function circle_points(F::CircleFun)
    circle_points(F.a,F.r,length(F.v))
end
