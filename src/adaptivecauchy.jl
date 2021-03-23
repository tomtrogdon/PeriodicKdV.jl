function Jacobi(a,b,n) # creates (n + 1) x (n + 1) Jacobi matrix
   SymTridiagonal([a(i) for i in 0:n],[b(i) for i in 0:n-1])
end

function dist(z,n) # check if inside Bernstein ellipse that tends to
    # [-1,1] as n -> ∞
    ρ = 1 + 1/n
    a = (ρ+1/ρ)/2.0
    b = (ρ-1/ρ)/2.0
    if real(z)^2/a^2 + imag(z)^2/b^2 <= 1
        return 1
    else
        return 0
    end
end


# This is all code to do adaptive QR with a tridiagonal matrix
import Base: complex

mutable struct QuadDiagonal  # length(diag) >= 3
   sub::Vector
   diag::Vector
   sup::Vector
   supsup::Vector
end

function Array(A::QuadDiagonal)
    b = (-1 => A.sub, 0 => A.diag, 1 => A.sup, 2 => A.supsup)
    diagm(b...)
end

function toBandedMatrix(A::QuadDiagonal)
    b = (-1 => A.sub, 0 => A.diag, 1 => A.sup, 2 => A.supsup)
    BandedMatrix(b,(length(A.diag),length(A.diag)))
end

function extend!(A::QuadDiagonal,v::Vector)
    append!(A.sub,v[1])
    append!(A.diag,v[2])
    append!(A.sup,v[3])
    append!(A.supsup,v[4])
end

function complex(A::QuadDiagonal)
    QuadDiagonal(A.sub |> complex, A.diag |> complex, A.sup |> complex, A.supsup |> complex) 
end

function complex!(A::QuadDiagonal)
    A.sub = A.sub |> complex
    A.diag = A.diag |> complex
    A.sup = A.sup |> complex
    A.supsup = A.supsup |> complex
end

function givens!(A::QuadDiagonal,j)
    if j > length(A.diag) -1
        @error "Matrix too small"
        return [1,0]
    end
    a = A.diag[j]
    b = A.sub[j]
    nom = sqrt(abs2(a)+abs2(b))
    c = a/nom; s = b/nom
    A.diag[j] = nom
    A.sub[j] = 0.0
    if j < n-1
        if A.supsup[j] != 0.0
        @warn "Structure issue, factorization not correct"
        end
        R = [conj(c) conj(s); -s c]*[A.sup[j] 0.0; A.diag[j+1] A.sup[j+1]]
        A.sup[j] = R[1,1]
        A.supsup[j] = R[1,2]
        A.diag[j+1] = R[2,1]
        A.sup[j+1] = R[2,2]
    else
        R = [conj(c) conj(s); -s c]*[A.sup[j] ; A.diag[j+1]]
        A.sup[j] = R[1]
        A.diag[j+1] = R[2]
    end
    [c,s]
end

function givens!(A::QuadDiagonal,v::Vector,j)
    n = length(A.diag)
    if j > n -1 || length(v) != n
        @error "Matrix wrong size"
        return A,v
    end
    a = A.diag[j]
    b = A.sub[j]
    nom = sqrt(abs2(a)+abs2(b))
    c = a/nom; s = b/nom
    A.diag[j] = nom
    A.sub[j] = 0.0
    if j < n-1
        if A.supsup[j] != 0.0
            @warn "Structure issue, factorization not correct"
        end
        R = [conj(c) conj(s); -s c]*[A.sup[j] 0.0 v[j]; A.diag[j+1] A.sup[j+1] v[j+1]]
        A.sup[j] = R[1,1]
        A.supsup[j] = R[1,2]
        A.diag[j+1] = R[2,1]
        A.sup[j+1] = R[2,2]
        v[j] = R[1,3]
        v[j+1] = R[2,3]
    else
        R = [conj(c) conj(s); -s c]*[A.sup[j] v[j]; A.diag[j+1] v[j+1]]
        A.sup[j] = R[1,1]
        A.diag[j+1] = R[2,1]
        v[j] = R[1,2]
        v[j+1] = R[2,2]
    end
    [c,s]
end

function QuadDiagonal(J::SymTridiagonal)
   QuadDiagonal(diag(J,-1) |> complex, diag(J) |> complex, diag(J,1) |> complex, fill(0.0im,length(diag(J))-2)) 
end

function pad(v,n)
    if length(v) == n
        return v
    elseif length(v) < n
        return vcat(v,zeros(typeof(v[1]),n-length(v)))
    else
        return v[1:n]
    end
end

function cauchy_off(a,b,nn,z)
    n = 3
    A = Jacobi(a,b,n) - z*I |> QuadDiagonal;
    v = fill(0.0im,n+1)
    v[1] = 1.0/(2im*pi)
    cs = [0im,0im]
    for i = 1:n
        cs = givens!(A,v,i)
    end
    n += 1
    while abs(v[end]) > 1e-16 && n < 10000
        c = cs[1]; s = cs[2]
        #R = [conj(c) conj(s); -s c]*[0im, b(n)]
        #R = b(n)*[conj(s), c]
        extend!(A,[b(n),a(n+1)-z, b(n)*c , b(n)*conj(s)])
        append!(v,0.0)
        cs = givens!(A,v,n)
        n += 1
    end
    pad(toBandedMatrix(A)\v,nn+1)
end

function cauchy(a,b,seed,n,z)
    if  dist(z,n) == 0   # the criterion for changing.
        # Non-adaptive method
        #m = 2; #over sampling, should be done adaptively
        #v = fill(0.0im,m*n+1)
        #v[1] = 1.0/(2im*pi)
        #c = ((Jacobi(a,b,m*n) - z*I)\v)
        # Adaptive method
        c = cauchy_off(a,b,n,z)
    else
        c = fill(0.0im,n+3)
        c[1] = seed(z);
        c[2] = z*c[1] - a(0)*c[1] + 1/(2im*pi)
        c[2] = c[2]/b(0)
        for j = 1:n-1 # compute c_n
            c[j+2] = z*c[j+1] - a(j)*c[j+1] - b(j-1)*c[j]
            c[j+2] /= b(j)
        end
    end
    c[1:n+1]
end

function cauchy(a,b,seed,n,z::Vector)
    vcat(map(zz -> cauchy(a,b,seed,n,zz) |> transpose, z)...)
end

### Still needed?

# function cauchyW(I1,n,I2,m)
#     a1 = I1[1]; b1 = I1[2]
#     M1 = x -> (b1-a1)/2*x .+ (b1+a1)/2
#     iM1 = x -> 2/(b1-a1)*x .- (b1+a1)/(b1-a1)
    
#     a2 = I2[1]; b2 = I2[2]
#     M2 = x -> (b2-a2)/2*x .+ (b2+a2)/2
#     #iM2 = x -> 2/(b2-a2)*x .- (b2+a2)/(b2-a2)
    
#     cauchy(aW,bW,seedW,n-1,iM1(M2(Ugrid(m))))
# end

# function cauchyW(I1,n,o::Integer,m)
#     a1 = I1[1]; b1 = I1[2]
#     M1 = x -> (b1-a1)/2*x .+ (b1+a1)/2
#     iM1 = x -> 2/(b1-a1)*x .- (b1+a1)/(b1-a1)
#     cauchy(aW,bW,seedW,n-1, Ugrid(m) .+ 1im*o*eps )
# end

# function cauchyV(I1,n,I2,m)
#     a1 = I1[1]; b1 = I1[2]
#     M1 = x -> (b1-a1)/2*x .+ (b1+a1)/2
#     iM1 = x -> 2/(b1-a1)*x .- (b1+a1)/(b1-a1)
    
#     a2 = I2[1]; b2 = I2[2]
#     M2 = x -> (b2-a2)/2*x .+ (b2+a2)/2
#     #iM2 = x -> 2/(b2-a2)*x .- (b2+a2)/(b2-a2)
    
#     cauchy(aV,bV,seedV,n-1,iM1(M2(Ugrid(m))))
# end

# function cauchyV(I1,n,o::Integer,m)
#     a1 = I1[1]; b1 = I1[2]
#     M1 = x -> (b1-a1)/2*x .+ (b1+a1)/2
#     iM1 = x -> 2/(b1-a1)*x .- (b1+a1)/(b1-a1)
#     cauchy(aV,bV,seedV,n-1, Ugrid(m) .+ 1im*o*eps )
# end

## Basic Cauchy
function Cauchy(WI1::WeightedInterval,WI2::WeightedInterval,n,m,o = 1)
    if WI1 == WI2
        a1 = WI1.a; b1 = WI1.b
        cauchy(WI1.W.a,WI1.W.b,WI1.W.seed,n-1, Ugrid(m) .+ 1im*o*eps )
    else
        a1 = WI1.a; b1 = WI1.b
        iM1 = iM(a1,b1)
    
        a2 = WI2.a; b2 = WI2.b
        M2 = M(a2,b2)
    
        return cauchy(WI1.W.a,WI1.W.b,WI1.W.seed,n-1,iM1(M2(Ugrid(m))))
    end
    
end

function Cauchy(C1::Cut,C2::Cut,n,m,o = 1)
    Cauchy(C1.W,C2.W,n,m,o)
end

function Cauchy(C1::Cut,WI2::WeightedInterval,n,m,o = 1)
    Cauchy(C1.W,WI2,n,m,o)
end
function Cauchy(WI1::WeightedInterval,C2::Cut,n,m,o = 1)
    Cauchy(WI1,C2.W,n,m,o)
end

### Basic CauchyChop
function CauchyChop(WI1::WeightedInterval,WI2::WeightedInterval,n,m,o,tol)
   chop(Cauchy(WI1,WI2,n,m,o) |> MatrixOperator, tol) 
end

function CauchyChop(C1::Cut,C2::Cut,n,m,o,tol)
   CauchyChop(C1.W,C2.W,n,m,o,tol)
end

function CauchyChop(WI1::WeightedInterval,C2::Cut,n,m,o,tol)
   CauchyChop(WI1,C2.W,n,m,o,tol)
end

function CauchyChop(C1::Cut,WI2::WeightedInterval,n,m,o,tol)
   CauchyChop(C1,WI2,n,m,o,tol)
end

### Cauchy on arrays of cuts or weighted intervals --- outputs a large matrix
function Cauchy(WI1::WeightedInterval,WI2::Array,n,m,o = 1)
   vcat(map(x -> Cauchy(WI1,x,n,m,o),WI2)...)
end

function Cauchy(C1::Cut,WI2::Array,n,m,o = 1)
   vcat(map(x -> Cauchy(C1.W,x,n,m,o),WI2)...)
end

function Cauchy(WI1::Array,WI2::Array,n,m,o = 1)
   hcat(map(x -> Cauchy(x,WI2,n,m,o),WI1)...)
end

### CauchyChop on arrays of cuts or weighted intervals --- outputs a BlockOperator

function CauchyChop(WI1::WeightedInterval,WI2::Array,n,m,o,tol)
   map(x -> CauchyChop(WI1,x,n,m,o,tol),WI2)
end

function CauchyChop(WI1::WeightedInterval,WI2::Array,n,m::Array,o,tol)
   map( (x,y) -> CauchyChop(WI1,x,n,y,o,tol),WI2,m)
end

function CauchyChop(C1::Cut,WI2::Array,n,m,o,tol)
   map(x -> CauchyChop(C1.W,x,n,m,o,tol),WI2)
end

function CauchyChop(C1::Cut,WI2::Array,n,m::Array,o,tol)
   map( (x,y) -> CauchyChop(C1.W,x,n,y,o,tol),WI2,m)
end

function CauchyChop(WI1::Array,WI2::Array,n,m,o,tol)
   hcat(map(x -> CauchyChop(x,WI2,n,m,o,tol),WI1)...) |> BlockOperator
end

function CauchyChop(WI1::Array,WI2::Array,n::Array,m::Array,o,tol)
   hcat(map( (x,y) -> CauchyChop(x,WI2,y,m,o,tol),WI1,n)...) |> BlockOperator
end

## TODO: Implement Cauchy(F::WeightFun,z::Complex{Float64},o) and Cauchy(F::Vector{WeightFun},z::Complex{Float64},o)

function Cauchy(F::WeightFun,z::Complex{Float64},o = 1)
   (cauchy(F.W.W.a,F.W.W.b,F.W.W.seed,length(F.cs)-1,[iM(F.W.a,F.W.b)(z+1im*o*eps)])*F.cs)[1]
end

function Cauchy(F::Vector{WeightFun},z::Complex{Float64},o = 1)
   map( x -> Cauchy(x,z,o), F) |> sum
end

function Cauchy(F,z::Float64,o = 1)
   Cauchy(F,complex(z),o)
end
