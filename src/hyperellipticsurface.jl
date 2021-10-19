mutable struct HyperellipticSurface
    #A::Array{Complex{Float64},2}  # Matrix of a-periods of basis differentials
    gaps::Array{Float64,2}  # gaps between cuts, suppose that these number are all positive, one band starts at zero
    g::Int64 # genus
    D::Array{Float64,2} # a simple divisor of points in the gaps, second entry \pm 1 for sheet
    Ω0::Vector{Complex{Float64}} # Abel map of D - D' where D' = gaps[:,1]
    Ωx::Vector{Complex{Float64}}
    Ωt::Vector{Complex{Float64}}
    E # a special constant
    α1
end

function HyperellipticSurface(q0::Function,L::Float64,n=100,tol=1e-14,m=200,trunctol=1e-14,invtol=.4,new=[true,true])
    gaps, zs, α1 = ScatteringData(q0,L,n,tol,trunctol)
    k = size(gaps)[1]
    S = HyperellipticSurface(gaps,zs,α1,m,new)
    RefineHyperellipticSurface!(q0,L,S,invtol,n,tol,k,m,new)
    return S
end

function HyperellipticSurface(gaps,zs,α1,m=50,new=[true,true])
    ep = 0 #TODO: Probably a good idea to eliminate the need for this.
    r = (x,y) -> (x |> Complex |> sqrt)/(y |> Complex |> sqrt)
    pr = (x,y) -> (x |> Complex |> sqrt)*(y |> Complex |> sqrt)
    #p = (z) -> -1im*prod(map(r,gaps[:,1] .- z, gaps[:,2] .- z))/sqrt(-z |> Complex)
    #P = (z) -> prod(map(pr, z .- gaps[:,1],  z .- gaps[:,2]))*1im*sqrt(-z |> Complex)
    bands = copy(gaps)
    bands[:,2] = gaps[:,1]
    bands[2:end,1] = gaps[1:end-1,2]
    bands[1,1] = 0.0

    # TODO: The biggest bottle neck in this whole thing is probably the evaluation of F.
    # Right now it is evaluated independently.  Should be able to iteration and use
    # F(j,k,z) to compute F(j',k',z) for neighboring j',k'
    function F(j,k,z)
        out = -1im/sqrt(-z |> Complex)
        for i = 1:size(gaps)[1]
            if i != j && i != k
                out *= r(z - gaps[i,1], z - gaps[i,2])
            end
        end
        if j == k
            return -out*pi*1im
        else
            return -1im*(gaps[k,2]-gaps[k,1])/2*out/pr(z - gaps[j,1], z - gaps[j,2])*pi
        end
    end


#     function G(j,k,z)
#         out = -1im/sqrt(gaps[end,end]-z |> Complex)
#         for i = 1:size(bands)[1]
#             if i != j && i != k
#                 out *= r(z - bands[i,2], z - bands[i,1])
#             end
#         end
#         if j == k
#             return -out*pi*1im
#         else
#             return 1im*(bands[k,2]-bands[k,1])/2*out/pr(z - bands[j,1], z - bands[j,2])*pi
#         end
#     end

    function γ(z)
        out = -1im/sqrt(-z |> Complex)
        for i = 1:size(gaps)[1]
            out *= r(z - gaps[i,1], z - gaps[i,2])
        end
        if imag(z) >= 0.0
            return out
        else
            return -out
        end
    end

	function γ(j,z)
        out = -1im/sqrt(-z |> Complex)
        for i = 1:size(gaps)[1]
			if i != j
            	out *= r(z - gaps[i,1], z - gaps[i,2])
			end
        end
        if imag(z) >= 0.0
            return out
        else
            return -out
        end
    end

    g = size(gaps)[1]
    A = zeros(Complex{Float64},g,g);

    if new[1]
		println("a-cycles: new method")
        as = (gaps[:,1] + gaps[:,2])/2 |> complex
        dist = (gaps[2:end,2] - gaps[1:end-1,1])/2 |> minimum # no need to be uniform here
	    dist = min(dist, gaps[1,1]/2)
        rads = (gaps[:,2] - gaps[:,1])/2 .+ dist
        cfuns = (a,rad) -> CircleFun(z -> map(γ,z), a, rad, m)
        γs = map(cfuns,as,rads)

        for i = 1:g
            b = gaps[i,1]
            divfun = ff -> CircleFun(ff.a,ff.r,ff.v./(circle_points(ff) .- b))
            df = map(divfun,γs)
            A[i,:] = map(Integrate,df)
        end
    else
        for i = 1:g
            for j = 1:g
                if i == j
                    A[i,i] = -2*(DefiniteIntegral(transformT,gaps[i,1],gaps[i,2],m)*(z -> F(i,i,z |> complex)))
                else
                    A[j,i] = -2*(DefiniteIntegral(transformV,gaps[i,1],gaps[i,2],m)*(z -> F(j,i,z |> complex)))
                end
            end
        end
    end

	#display(A)

#     tB = zeros(Complex{Float64},g,g);
#     for i = 1:g
#         for j = 1:g
#             if i == j
#                 tB[i,i] = 2*(DefiniteIntegral(transformT,bands[i,1],bands[i,2],m)*(z -> G(i,i,z+1im*ep)))
#             else
#                 tB[j,i] = 2*(DefiniteIntegral(transformW,bands[i,1],bands[i,2],m)*(z -> G(j,i,z+1im*ep)))
#             end
#         end
#     end

#     B = copy(tB)
#     B[:,1] = tB[:,1]
#     for j = 2:g
#        B[:,j] = B[:,j] + B[:,j-1]
#     end
#     B = 2im*pi*(A\B)  #Riemann Matrix

    Ωx = zeros(Complex{Float64},g);
	if new[1]
		for i = 1:g
		    b = gaps[i,1]
		    divfun = ff -> CircleFun(ff.a,ff.r,(ff.v./(circle_points(ff) .- b)).*sqrt.(circle_points(ff)))
		    df = map(divfun,γs)
		    Ωx[i] = -1im*sum(map(Integrate,df))
		end
	else
		for i = 1:g
        	for j = 1:g
            	if i == j
                	Ωx[i] += 2im*(DefiniteIntegral(transformT,gaps[i,1],gaps[i,2],m)*(z -> sqrt(z)*F(i,i,z+1im*ep)))
            	else
                	Ωx[j] += 2im*(DefiniteIntegral(transformV,gaps[i,1],gaps[i,2],m)*(z -> sqrt(z)*F(j,i,z+1im*ep)))
            	end
        	end
    	end
	end
    Ωx = -2*(A\Ωx)

    # TODO: For efficiency E should be a scalar, but this serves as a check,
    # E == E[1]*fill(1.0,g)

    E = zeros(Complex{Float64},g);
	if new[1]
		for i = 1:g
		    b = gaps[i,1]
		    divfun = ff -> CircleFun(ff.a,ff.r,(ff.v./(PeriodicKdV.circle_points(ff) .- b)).*PeriodicKdV.circle_points(ff))
		    df = map(divfun,γs)
		    E[i] = -0.5*sum(Ωx.*map(PeriodicKdV.Integrate,df))
		end
	else
    	for i = 1:g
        	for j = 1:g
            	if i == j
                	E[j] += Ωx[i]*(DefiniteIntegral(transformT,gaps[i,1],gaps[i,2],m)*(z -> z*F(i,i,z+1im*ep)))
            	else
                	E[j] += Ωx[i]*(DefiniteIntegral(transformV,gaps[i,1],gaps[i,2],m)*(z -> z*F(j,i,z+1im*ep)))
            	end
        	end
    	end
	end


    Ωt = zeros(Complex{Float64},g);
	if new[1]
		for i = 1:g
			b = gaps[i,1]
			divfun = ff -> CircleFun(ff.a,ff.r,(ff.v./(circle_points(ff) .- b)).*sqrt.(circle_points(ff)).^3)
			df = map(divfun,γs)
			Ωt[i] = -4im*sum(map(Integrate,df))
		end
	else
		for i = 1:g
        	for j = 1:g
            	if i == j
                	Ωt[i] += 8im*(DefiniteIntegral(transformT,gaps[i,1],gaps[i,2],m)*(z -> z^(3/2)*F(i,i,z+1im*ep)))
            	else
                	Ωt[j] += 8im*(DefiniteIntegral(transformV,gaps[i,1],gaps[i,2],m)*(z -> z^(3/2)*F(j,i,z+1im*ep)))
            	end
        	end
    	end
	end
    E -= Ωt/4
    E /= 2im*pi
    Ωt = -2*(A\Ωt)

	function J(n,λ)
        if λ > 1 && n == 0
            return 1.0
        elseif λ > 1
            return 0
        end

        if λ < -1
            return 0.0
        end

        if n == 0
            1 - acos(λ)/pi
        else
            - sqrt(2)/pi*sin(n*acos(λ))/n
        end
    end

    function J(f,a,b,n,λ)
        cs = transformT(map(f,M(a,b)(Ugrid(n))))
        #display(cs[end] |> abs)
        out = 0im
        for i = 1:length(cs)
            out += cs[i]*J(i-1,iM(a,b)(λ))
        end
        cs[1] - out
    end

	function J(v,a,b,λ)
	    cs = transformT(v)
	    #display(cs[end] |> abs)
	    out = 0im
	    for i = 1:length(cs)
	        out += cs[i]*J(i-1,iM(a,b)(λ))
	    end
	    cs[1] - out
	end


	if new[2]
		println("Abel: new method")
		gen_abel_pts = (a,b,λ) -> M(a,b)(Ugrid(m))
		abel_pts = map(gen_abel_pts,gaps[:,1],gaps[:,2],zs[:,1])
		abel_γ =  copy(abel_pts)
		for i = 1:g
			abel_γ[i] = map(z -> γ(i,z), abel_pts[i])
		end

		function Abelvec(n,j,k,λ) # integrate differential k over part of gap j
		    a = gaps[j,1]
		    b = gaps[j,2]
		    if k == j
		        -J(xx -> F(j,k,xx), a, b, n, λ)
		    else
		        vals = copy(abel_γ[j])
		        #vals .*= map(r,abel_pts[j] .- b, abel_pts[j] .- a)
		        vals .*= (abel_pts[j] .- a)./(abel_pts[j] .- gaps[k,1])
		        1im*pi*J(vals,a,b,λ)
		    end
		end

		function Abelvec(n,k,λ) # integrate differential k over part of gap j
		    j = 1
		    for jj = 1:g
		        if gaps[j,1] <= λ <= gaps[j,2]
		            break
		        end
		        j += 1
		    end
		    if j == g + 1
		        #@warn "Not in a gap"
		        return 0
		    else
		        return Abelvec(n,j,k,λ)
		    end
		end
	else
    	function Abel(n,j,k,λ) # integrate differential k over part of gap j
        	a = gaps[j,1]
        	b = gaps[j,2]
        	if k == j
            	-J(xx -> F(j,k,xx), a, b, n, λ)
        	else
            	#println(2/(gaps[j,2]-gaps[j,1]))
            	-2/(gaps[j,2]-gaps[j,1])*J( xx -> F(k,j,xx)*(xx - gaps[j,1]), a, b, n, λ)
        	end
    	end

    	function Abel(n,k,λ) # integrate differential k over part of gap j
        	j = 1
        	for jj = 1:g
            	if gaps[j,1] <= λ <= gaps[j,2]
                	break
            	end
            	j += 1
        	end
        	if j == g + 1
            	#@warn "Not in a gap"
            	return 0
        	else
            	return Abel(n,j,k,λ)
        	end
    	end
		println("Abel: old method")
	end

	if new[2]
		abel =  map( k -> -zs[1,2]*Abelvec(m,k,zs[1,1]), 1:g)
		for j = 2:g
			 abel += map( k -> -zs[j,2]*Abelvec(m,k,zs[j,1]), 1:g)
		end
	else
		abel =  map( k -> -zs[1,2]*Abel(m,k,zs[1,1]), 1:g)
    	for j = 2:g
        	abel += map( k -> -zs[j,2]*Abel(m,k,zs[j,1]), 1:g)
    	end
	end
    #display((A,abel))
    abel = (A\abel)*(2im*pi)

    HyperellipticSurface(gaps,g,zs,abel,Ωx,Ωt,E,α1)

end

struct BakerAkhiezerFunction
    WIm::Vector{WeightedInterval}
    WIp::Vector{WeightedInterval}
    Ω::Function
    E::Complex{Float64}
    α1
    Cp # Compressed Cauchy matrix
    Cm
    ns
    tol
    iter
end

Skew = θ -> [0 exp(θ); exp(-θ) 0]
Skewx = (θ,θx) -> [0 θx*exp(θ); -θx*exp(-θ) 0]
gp = (x,θ) -> Cut(Skew(θ), x)
gm = (x,θ) -> Cut(Skew(-θ), x)

## Swap definitions?
gpx = (x,θ,θx) -> Cut(Skewx(θ,θx), x)
gmx = (x,θ,θx) -> Cut(Skewx(-θ,-θx), x)


function bernsteinρ(a,b,r)
    rr = 2*r/(b-a) + 1.0
    ρ =  rr - sqrt(rr^2 - 1)
end

function bernsteinρ(v::Array)
    u = copy(v[:,1])
    rad = min(v[2,1]-v[1,2],2*v[1,1])
    u[1] = bernsteinρ(v[1,1],v[1,2],rad/2)
    for i = 2:length(u)-1
        rad = min(v[i,1]-v[i-1,2],v[i+1,1]-v[i,2])
        u[i] = bernsteinρ(v[i,1],v[i,2],rad/2)
    end
    rad = v[end,1]-v[end-1,2]
    u[end] = bernsteinρ(v[end,1],v[end,2],rad/2)
    u
end

fff = (ϵ,ρ,c,k) -> max(ρ < 1e-16 ? 0 : convert(Int64, (log(1/ϵ) + log(c/(1-ρ)))/log(1/ρ) |> ceil),k)

function choose_order(gaps::Array,ϵ,c,k)
    map(z -> fff(ϵ,z,c,k), bernsteinρ(gaps)) .+ 2
end

function BakerAkhiezerFunction(S::HyperellipticSurface,n::Int64,tols = [2*1e-14,false],iter = 100)
    zgaps_neg = hcat(- sqrt.(S.gaps[:,2]) |> reverse, - sqrt.(S.gaps[:,1]) |> reverse)
    zgaps_pos = hcat( sqrt.(S.gaps[:,1]) , sqrt.(S.gaps[:,2]) )
    #zzs_pos = sqrt.(zs)
    #zzs_neg = -sqrt.(zs) |> reverse;
    fV = (x,y) -> WeightedInterval(x,y,chebV)
    fW = (x,y) -> WeightedInterval(x,y,chebW)
    Ω = (x,t) -> -S.Ωx*x - S.Ωt*t - S.Ω0

    WIm = map(fW,zgaps_neg[:,1],zgaps_neg[:,2])
    WIp = map(fV,zgaps_pos[:,1],zgaps_pos[:,2])

    Ωs = Ω(0.0,0.0)
    RHP = vcat(map(gm,WIm,Ωs |> reverse),map(gp,WIp,Ωs));
    ns = fill(n,WIp |> length)

#     #lens = abs.(zgaps[:,1] - zgaps[:,2])


#     f = x -> convert(Int,ceil(10 + 10/x^2))
#     ns = map(f,zgaps_pos[:,1])
    ns = vcat(ns |> reverse, ns)

    CpBO = CauchyChop(RHP,RHP,ns,ns,1,tols[2])
    CmBO = CauchyChop(RHP,RHP,ns,ns,-1,tols[2])


    #println("Effective rank of Cauchy operator = ",effectiverank(CpBO))
    #println("Maximum rank of Cauchy operator = ", (2*S.g)^2*n )


    return BakerAkhiezerFunction(WIm,WIp,Ω,S.E[1],S.α1,CpBO,CmBO,ns,tols[1],iter)
end

function BakerAkhiezerFunction(S::HyperellipticSurface,c::Float64,tols = [2*1e-14,false],iter = 100,K=0,show_flag=false)
    zgaps_neg = hcat(- sqrt.(S.gaps[:,2]) |> reverse, - sqrt.(S.gaps[:,1]) |> reverse)
    zgaps_pos = hcat( sqrt.(S.gaps[:,1]) , sqrt.(S.gaps[:,2]) )
    #zzs_pos = sqrt.(zs)
    #zzs_neg = -sqrt.(zs) |> reverse;
    fV = (x,y) -> WeightedInterval(x,y,chebV)
    fW = (x,y) -> WeightedInterval(x,y,chebW)
    Ω = (x,t) -> -S.Ωx*x - S.Ωt*t - S.Ω0

    WIm = map(fW,zgaps_neg[:,1],zgaps_neg[:,2])
    WIp = map(fV,zgaps_pos[:,1],zgaps_pos[:,2])

    Ωs = Ω(0.0,0.0)
    RHP = vcat(map(gm,WIm,Ωs |> reverse),map(gp,WIp,Ωs));
    ns = choose_order(zgaps_pos,tols[1],c,K)
    if show_flag
    	println(ns)
    end


#     #lens = abs.(zgaps[:,1] - zgaps[:,2])


#     f = x -> convert(Int,ceil(10 + 10/x^2))
#     ns = map(f,zgaps_pos[:,1])
    ns = vcat(ns |> reverse, ns)

    CpBO = CauchyChop(RHP,RHP,ns,ns,1,tols[2])
    CmBO = CauchyChop(RHP,RHP,ns,ns,-1,tols[2])


    #println("Effective rank of Cauchy operator = ",effectiverank(CpBO))
    #println("Maximum rank of Cauchy operator = ", (2*S.g)^2*n )


    return BakerAkhiezerFunction(WIm,WIp,Ω,S.E[1],S.α1,CpBO,CmBO,ns,tols[1],iter)
end

function BakerAkhiezerFunction(S::HyperellipticSurface,c::Array,tol = 2*1e-14,iter = 100)
    zgaps_neg = hcat(- sqrt.(S.gaps[:,2]) |> reverse, - sqrt.(S.gaps[:,1]) |> reverse)
    zgaps_pos = hcat( sqrt.(S.gaps[:,1]) , sqrt.(S.gaps[:,2]) )
    #zzs_pos = sqrt.(zs)
    #zzs_neg = -sqrt.(zs) |> reverse;
    fV = (x,y) -> WeightedInterval(x,y,chebV)
    fW = (x,y) -> WeightedInterval(x,y,chebW)
    Ω = (x,t) -> -S.Ωx*x - S.Ωt*t - S.Ω0

    WIm = map(fW,zgaps_neg[:,1],zgaps_neg[:,2])
    WIp = map(fV,zgaps_pos[:,1],zgaps_pos[:,2])

    Ωs = Ω(0.0,0.0)
    RHP = vcat(map(gm,WIm,Ωs |> reverse),map(gp,WIp,Ωs));
    ns = c # choose_order(zgaps_pos,tol,c)
#     #lens = abs.(zgaps[:,1] - zgaps[:,2])
#     f = x -> convert(Int,ceil(10 + 10/x^2))
#     ns = map(f,zgaps_pos[:,1])
    ns = vcat(ns |> reverse, ns)
    CpBO = CauchyChop(RHP,RHP,ns,ns,1,tol)
    CmBO = CauchyChop(RHP,RHP,ns,ns,-1,tol)
    #println("Effective rank of Cauchy operator = ",effectiverank(CpBO))
    #println("Maximum rank of Cauchy operator = ", (2*S.g)^2*n )
    return BakerAkhiezerFunction(WIm,WIp,Ω,S.E[1],S.α1,CpBO,CmBO,ns,tol,iter)
end

function (BA::BakerAkhiezerFunction)(x,t,tol = BA.tol)
    ns = BA.ns
    Ωs = BA.Ω(x,t)
    Ωsx = BA.Ω(1.0,0) - BA.Ω(0.0,0)

    RHP = vcat(map(gm,BA.WIm,Ωs |> reverse),map(gp,BA.WIp,Ωs));
    RHPx = vcat(map(gmx,BA.WIm, Ωs |> reverse, Ωsx |> reverse),map(gpx,BA.WIp,Ωs,Ωsx));

    m = length(RHP)
    #Z = ZeroOperator(n)
    #ZZ = BlockOperator(fill(Z,m,m))
    ZZ = BlockZeroOperator(ns)
    C⁺ = vcat(hcat(BA.Cp,ZZ),hcat(ZZ,BA.Cp))


    D11 = DiagonalBlockOperator( [-RHP[j].J[1,1] for j = 1:m], ns)
    D21 = DiagonalBlockOperator( [-RHP[j].J[1,2] for j = 1:m], ns)
    D12 = DiagonalBlockOperator( [-RHP[j].J[2,1] for j = 1:m], ns)
    D22 = DiagonalBlockOperator( [-RHP[j].J[2,2] for j = 1:m], ns)

    D11x = DiagonalBlockOperator( [-RHPx[j].J[1,1] for j = 1:m], ns)
    D21x = DiagonalBlockOperator( [-RHPx[j].J[1,2] for j = 1:m], ns)
    D12x = DiagonalBlockOperator( [-RHPx[j].J[2,1] for j = 1:m], ns)
    D22x = DiagonalBlockOperator( [-RHPx[j].J[2,2] for j = 1:m], ns)

    p = [1,m+1]
    for i = 2:m
        append!(p,i)
        append!(p,m+ i)
    end

    JC⁻ = vcat(hcat(D11*BA.Cm,D12*BA.Cm),hcat(D21*BA.Cm,D22*BA.Cm))
    JxC⁻ = vcat(hcat(D11x*BA.Cm,D12x*BA.Cm),hcat(D21x*BA.Cm,D22x*BA.Cm))
    S = C⁺ + JC⁻
    # notes for x derivative of solution
    # S*u = b
    # Sx*u + S*ux = bx
    # S*ux = bx - Sx*u,   Sx = JxC-

    Sp = permute(S,p)
    JxCp = permute(JxC⁻,p)

    D = TakeDiagonalBlocks(Sp,2)

    #### fix this.
    ind = vcat(ns,ns)
    dim = sum(ind)
    ind = ind |> gapstoindex
    fine_ind = indextogaps(ind)[p] |> gapstoindex
    coarse_ind = fine_ind |> indextogaps
    coarse_ind = coarse_ind[1:2:end] + coarse_ind[2:2:end] |> gapstoindex

    PrS = x -> BlockVector(fine_ind,D\BlockVector(coarse_ind,x))

    b = BlockVector(ind,fill(0.0im,2*dim))
    bx =  BlockVector(ind,fill(0.0im,2*dim))

    #### Do this using block vectors

    for j in 1:m
        b[j] = b[j] .+ (RHP[j].J[1,1] + RHP[j].J[2,1] - 1)
        b[j+m] = b[j+m] .+ (RHP[j].J[2,2] + RHP[j].J[1,2] - 1)
        bx[j] = bx[j] .+ (RHPx[j].J[1,1] + RHPx[j].J[2,1])
        bx[j+m] = bx[j+m] .+ (RHPx[j].J[2,2] + RHPx[j].J[1,2])

#         rhs[(j-1)*n+1:j*n] = rhs[(j-1)*n+1:j*n] .+ (RHP[j].J[1,1] + RHP[j].J[2,1] - 1)
#         rhs[m*n+(j-1)*n+1:m*n+j*n] = rhs[m*n+(j-1)*n+1:m*n+j*n] .+ (RHP[j].J[2,2] + RHP[j].J[1,2] - 1)
#         rhsx[(j-1)*n+1:j*n] = rhsx[(j-1)*n+1:j*n] .+ (RHPx[j].J[1,1] + RHPx[j].J[2,1])
#         rhsx[m*n+(j-1)*n+1:m*n+j*n] = rhsx[m*n+(j-1)*n+1:m*n+j*n] .+ (RHPx[j].J[2,2] + RHPx[j].J[1,2])
    end

    #b = BlockVector(rhs,n)
    #bx = BlockVector(rhsx,n)
    bp = permute(b,p)
    bpx = permute(bx,p)

    Op = x -> PrS(Sp*x)
    #return (D,coarse_ind,fine_ind,bp)
    out = GMRES_quiet(Op,PrS(bp),⋅, tol,BA.iter)
    solp = out[2][1]*out[1][1]
    for j = 2:length(out[2])
        solp += out[2][j]*out[1][j]
    end

    outx = GMRES_quiet(Op,PrS(bpx-JxCp*solp),⋅, tol,BA.iter)
    solpx = outx[2][1]*outx[1][1]
    for j = 2:length(outx[2])
        solpx += outx[2][j]*outx[1][j]
    end

    #Old way
    #solp = +([out[2][j]*out[1][j] for j = 1:length(out[2])]...)

    sol = ipermute(solp,p)
    solx = ipermute(solpx,p)

    f1 = [WeightFun(sol[j],RHP[j].W)  for j = 1:m ]
    f2 = [WeightFun(sol[m + j],RHP[j].W)  for j = 1:m ]
    f1x = [WeightFun(solx[j],RHP[j].W)  for j = 1:m ]
    f2x = [WeightFun(solx[m + j],RHP[j].W)  for j = 1:m ]

    Φ = (z,o) -> [Cauchy(f1,z,o) + 1.0, Cauchy(f2,z,o) + 1.0] |> transpose
    (Φ,f1,f2,f1x,f2x,BA.E)
    #(Φ,f1,f2)
end

function KdV(BA::BakerAkhiezerFunction,x,t)
    out = BA(x+6*BA.α1*t, t);
    2im*(-1/(2im*pi)*sum(map(x -> DomainIntegrateVW(x), out[4])) + BA.E) - BA.α1
    # I think this is right but I cannot justify it, +/- sign issue
    # It must have to do with the jumps and which sheet, etc.
end

function KdV(BA::BakerAkhiezerFunction,x,t,tol)
    out = BA(x+6*BA.α1*t, t, tol);
    2im*(-1/(2im*pi)*sum(map(x -> DomainIntegrateVW(x), out[4])) + BA.E) - BA.α1
    # I think this is right but I cannot justify it, +/- sign issue
    # It must have to do with the jumps and which sheet, etc.
end
