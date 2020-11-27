### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 23dab282-10cc-11eb-0ade-c7f191741071
using Random, Distributions, Plots, StatsBase, PlutoUI

# ╔═╡ 020b1160-22ef-11eb-2974-e9d551ada921
Random.seed!(123)

# ╔═╡ b93e3ee0-10c8-11eb-3db9-0f2cfebf6f7a
md"""
# Migration over a linear habitat

This notebook contains an extension of single deme population genetics to a linear habitat with migration.

"""

# ╔═╡ 6e53a382-10cc-11eb-0873-adb1d891562b
# please document! also don't rely on global variables! these should all be in one of the structs you define (either as fields, or lengths of arrays etc.). Unless they are input arguments to the main function call.
begin
		
	N = 50 #100     # deme level
	
	L = 250  # agent level
	
	p = 0.5  # initialization?
	
	t = 5000  # simulation 
	
	Vs = 1/2 #2  # deme level
	
	rm = 1.06 #1.025  # deme level
	
	α = 0.1 #0.25  # individual level
	
	θ = 12.5  # deme level
	
	K = 50 #125  # deme level
	
	μ = 0.000001  # deme level
	
	σ = sqrt(1/2)  # habitat level
	
	h = 0.1  # ?
	
	b = 0.1  # habitat level
	
	Dm = 250
end

# ╔═╡ 1eb68d20-17e2-11eb-3166-d7bff64901ed
# Could you document this a bit more? Why do we use these relations? I guess it's from Polechova & Barton? Why did they choose these kind of relations? This is the kind of stuff you will look at later and have forgotten why this was so, so some comments would be good.
begin
	num_demes(σ, Vs, α) = round(Int64, 10*4*sqrt(σ^2 * Vs)/α)
	num_loci(nd, b, α) = round(Int64, nd*b/α)
	num_ind(K, b, σ, h, Vs, rm) = round(Int64, K*((1.0 - b*σ)/(2*h^2*sqrt(Vs)*rm)))
end

# ╔═╡ 29c22850-17e2-11eb-06ed-492978d8d4ea
nd = num_demes(σ, Vs, α)

# ╔═╡ 3fe1c320-17e2-11eb-0759-d508c9753c04
nl = num_loci(nd, b, α)

# ╔═╡ 4093d510-17e2-11eb-2157-c5a8d98967c9
ni = num_ind(K, b, σ, h, Vs, rm)

# ╔═╡ 55fe6b60-10c9-11eb-18a6-e9963284f863
struct Agent{T}
    loci::Vector{T}
end

# ╔═╡ 48b0da00-10ca-11eb-2684-792d22df5a69
struct Deme{A,T}
    agents::Vector{A}
    K::Int64   
    θ::T      
end

# ╔═╡ 6d3d5040-1658-11eb-0e06-095bb2d6b3a5
struct Habitat{D}
    demes::Vector{D}
end

# ╔═╡ 17c43700-10cc-11eb-277b-bf13ce5f3505
begin
	randagent(d::Distribution, n) = Agent(rand(d, n))
	randagent(d::Distribution, n, N) = [Agent(rand(d, n)) for i=1:N]
end

# ╔═╡ 9920bd70-10ca-11eb-343a-4b068e4a4f16
begin
	Base.length(d::Deme) = length(d.agents)
	Base.length(a::Agent) = length(a.loci)
	Base.length(h::Habitat) = length(h.demes)
	Base.show(io::IO, d::Deme{A,T}) where {A,T} = 
		write(io, "Deme{$A,$T}(N=$(length(d)))")
	Base.show(io::IO, h::Habitat) = 
		write(io,"Habitat (n=$(length(h))):\n  $(join(string.(h.demes), "\n  "))")
	Base.push!(d::Deme, a::Agent) = push!(d.agents, a)
	Base.push!(a::Agent{T}, x::T) where T = push!(a.loci, x)
	Base.push!(h::Habitat, d::Deme) = push!(h.demes, d)
	Base.getindex(a::Agent, i) = a.loci[i]
	Base.getindex(d::Deme, i) = d.agents[i]
	Base.getindex(d::Deme, i, j) = d.agents[i][j]
	Base.rand(d::Deme) = rand(d.agents)
	Base.rand(d::Deme, n::Integer) = rand(d.agents, n)
	Base.sum(a::Agent) = sum(a.loci) 
end

# ╔═╡ e8008b00-10ca-11eb-135b-35d961573890
begin
	phenotype(agent::Agent) = sum(agent)
	phenotype(agent::Agent, zmean) = sum(agent) - zmean
end

# ╔═╡ ccb28910-12e7-11eb-0869-8be87d5a4ae6
function trait_mean(d::Deme)
	z = Float64[]
	for agent in d.agents
		push!(z,sum(agent))
	end
	sum(z)/length(d)
end		

# ╔═╡ 070648f2-10cb-11eb-2f40-3573a65cb627
begin
	emptycopy(a::Agent{L}) where L = Agent(L[])
	emptycopy(d::Deme{A}) where A = Deme(A[], d.K, d.θ)
	emptycopy(h::Habitat) = Habitat(emptycopy.(h.demes))
end

# ╔═╡ 7a674f20-3082-11eb-2956-01e2f8acb4be


# ╔═╡ 30074380-10cb-11eb-3600-f58f98f07fcd
function mate(a::Agent, b::Agent)
    newloci = similar(a.loci)
    for i in 1:length(a)
		newloci[i] = rand() < 0.5 ? a[i] : b[i] 
    end
	Agent(newloci)
end 

# ╔═╡ 3be1ee80-10cb-11eb-275b-b73888bb0061
function random_mating(d::Deme)
    newdeme = similar(d.agents)
    for i=1:length(d)
		newdeme[i] = mate(rand(d, 2)...)
	end
	Deme(newdeme, d.K, d.θ)
end 

# ╔═╡ f2877110-1223-11eb-1140-afde4ce678cd
function allelef(d::Deme)
	f = j->mapreduce(i->d[i,j], +, 1:length(d))/length(d)
	frequencies = map(f, 1:length(d[1]))
end

# ╔═╡ 74004a90-124d-11eb-1f21-39c3d9aac277
function allelefreq(d::Deme)
	freq = Vector{Float64}(undef,length(d[1]))
	for j in 1:length(d[1])
		s = 0
    	for i in d.agents
        	if i.loci[j] == α   # this is a global variable! avoid!
            	s += 1
			end
        end
		f = s/length(d)
		freq[j] = f
	end
	freq
end	

# ╔═╡ e8641300-1223-11eb-1d39-091974fefa58
heterozygosities(d::Deme, freqs=allelefreq(d)) = map(p->p*(1-p), freqs)

# ╔═╡ f71a1890-1223-11eb-1b02-4980890c068d
expected_heterozygosity(H₀, t, N) = ((1.0-1.0/N)^t)*H₀

# ╔═╡ b1739730-1223-11eb-28ff-f56636bc7119
function neutral_evolving_deme(d::Deme, ngen; fun=heterozygosities, trait_mean= trait_mean, allelefreq = allelefreq)
	stats = [fun(d)]
	tm = [trait_mean(d)]
	af = [allelefreq(d)]
	for n=1:ngen
		d = random_mating(d)
		push!(stats, fun(d))
		push!(tm, trait_mean(d))
		push!(af,allelefreq(d))
	end
	(stats=stats, tm = tm, af = af, deme=d, ngen=ngen)
end

# ╔═╡ f2213480-10cb-11eb-3243-59823606babf
d = Deme(randagent(DiscreteNonParametric([0,α],[p,1-p]), L, N), K, θ)

# ╔═╡ 6ff28710-12e8-11eb-2708-c90089167cd6
trait_mean(d)

# ╔═╡ 37ce3f12-124c-11eb-1cc0-05e92eaec860
allelef(d)

# ╔═╡ 1f5a1afe-124e-11eb-0819-b7a13a875f45
allelefreq(d)

# ╔═╡ 70c6d500-124e-11eb-3da9-e1bd16cc0a74
allelef(d) == allelefreq(d)

# ╔═╡ e76d7e00-1223-11eb-05c8-d71557c58345
heterozygosities(d)

# ╔═╡ b3c37b90-12dc-11eb-1628-6d4571848222
d[1]

# ╔═╡ b5933190-1223-11eb-091e-9dd0af287609
sim = neutral_evolving_deme(d, t)

# ╔═╡ c095cc12-1223-11eb-3927-019755ca8d18
begin
	Hₒ = map(mean, sim.stats)
	plot(Hₒ, grid=false, color=:black, label="\$H_o(t)\$")
	plot!(0:sim.ngen+1, 
		t->expected_heterozygosity(p*(1-p), t, N),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 839c8dfe-12e8-11eb-0a77-199fe2641df1
begin
	traitmean = map(mean, sim.tm)
	plot(traitmean, grid=false, color=:black, label="Random mating")
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 19191fe0-2abb-11eb-1f81-d591e973be82
begin
	af = sim.af
	scatter(af[length(af)], grid=false, color=:black, label="Random mating")
	xlabel!("Locus")
	ylabel!("Allele frequency after t generations")
end

# ╔═╡ 8f9ee77e-1233-11eb-25f9-df1d20c7e554
md"""
## Simulation including population regulation and stabilizing selection
"""

# ╔═╡ d07e5d80-10cb-11eb-2e72-3b1acc3996a1
md"""

The Malthusian fitness of an individual: $r(N,z) = r_e(N) + r_g(z)$, with 

$r_e(N) = r_m * (1 - N/K)$

$r_g(z) = -(z-θ)^2 * 1/2Vs$

"""

# ╔═╡ 41e53d00-10cb-11eb-3153-c1a586773c80
function malthusian_fitness(d::Deme,a::Agent,rm,Vs)
    N = length(d)
    z = sum(a)
    return rm*(1-(N/d.K))-((z-d.θ)^2)/(2*Vs)
end 

# ╔═╡ 18cf418e-10d0-11eb-36b1-1d5426742a23
function malthusian_fitness(d::Deme,rm,Vs)
    N = length(d)
    fitnesses = Float64[]
	for agent in d.agents
		z = sum(agent)
    	f = rm*(1-(N/d.K))-((z-d.θ)^2)/(2*Vs)
		push!(fitnesses, f)
	end
	fitnesses
end 

# ╔═╡ 0451e510-10d0-11eb-2d78-25fdaa8e5f5e
malthusian_fitness(d,rm,Vs)

# ╔═╡ 93b4d780-10cb-11eb-06f1-556745d3c323
md"""
Every generation, each individual $i$ produces a
number of offspring according to a Poisson distribution with mean the individual’s fitness $Exp(r_i)$.

"""



# ╔═╡ 8e06fc00-10cb-11eb-22d6-753424c41d8c
number_of_offspring(d::Deme,a::Agent,rm,Vs) = rand(Poisson(exp(malthusian_fitness(d::Deme,a::Agent,rm,Vs))))

# ╔═╡ bf49d880-10cd-11eb-28b4-4111e7c69a54
begin
	nof = []
	for agent in d.agents
	    push!(nof, number_of_offspring(d,agent,rm,Vs))
	end
	nof
end

# ╔═╡ 24b4a580-10d1-11eb-1bcb-d775fb8b26b9
function replicate(d::Deme{A}) where {A}
    newdeme = emptycopy(d)
    for agent in d.agents
        numb = number_of_offspring(d,agent,rm,Vs)
        for i in 1:numb
            push!(newdeme,agent)
        end
    end
    newdeme
end

# ╔═╡ a8fd1640-1259-11eb-259e-fdbd74535c88
length(d)

# ╔═╡ a13eaa3e-1259-11eb-2025-11cb34cf71f5
length(replicate(d))

# ╔═╡ cab06690-12de-11eb-07ae-232de34b05c3
md"""

##### Mating probability proportional to fitness

"""

# ╔═╡ 5017abf0-10d1-11eb-0f13-41564d34d322
md"""

Mates are sampled with replacement proportional to their fitness. The fitness $F_i$ of an individual i equals $Exp(r_i)$, where $r_i$ is the malthusian fitness as described earlier.

"""

# ╔═╡ 521567b0-1245-11eb-1c3a-e31e283d59c7
fitnesses = exp.(malthusian_fitness(d,rm,Vs))

# ╔═╡ 873f5780-1249-11eb-105c-d16a6160ec4f
sample(d.agents, weights(fitnesses),2)

# ╔═╡ 66e5d820-1243-11eb-2162-755e424f3c4c
function mating_extended(d::Deme{A}, fitnesses) where A
    newdeme = Vector{A}(undef, length(d))
    for i=1:length(d)
		newdeme[i] = mate(sample(d.agents, weights(fitnesses),2)...)
	end
	Deme(newdeme, d.K, d.θ)
end 

# ╔═╡ 447b6040-2079-11eb-0d3c-f7c17caae8dc
"""Mating function that I use now and is the most likely to be in accordance with how PnB implemented it"""
function mating_PnB(d::Deme{A}, rm, Vs) where A
	newdeme = Vector{A}(undef,0)
	fitnesses = exp.(malthusian_fitness(d,rm,Vs))
	for i=1:length(d)
		B1 = d.agents[i]
		noff = number_of_offspring(d,B1,rm,Vs)
		B2 = sample(d.agents, weights(fitnesses))
		#B2 = rand(d.agents)
		#child = mate(B1,B2)
		for c in 1:noff 
				push!(newdeme, mate(B1,B2))
		end
	end
	Deme(newdeme, d.K, d.θ)
end 		

# ╔═╡ 1ec87310-1232-11eb-0bbc-895d245932fc
function evolving_deme_extended(d::Deme, ngen, fitnesses; fun=heterozygosities)
	stats = [fun(d)]
	for n=1:ngen
		#d = mating_extended(d, fitnesses)
		d = mating_PnB(d, rm, Vs)
		push!(stats, fun(d))
	end
	(stats=stats, deme=d, ngen=ngen)
end

# ╔═╡ bbe52202-1233-11eb-06d4-7f58dcfed50e
sim_extended = evolving_deme_extended(d, t, fitnesses)

# ╔═╡ c6b23d80-1247-11eb-2ef5-a7295f4aff4f
begin
	Hₒ_ext = map(mean, sim_extended.stats)
	plot(Hₒ_ext, color=:black, label="\$H_o(t)\$")
	plot!(0:sim_extended.ngen+1, 
		t->expected_heterozygosity(p*(1-p), t, N),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 052f6f12-12d9-11eb-2c06-5fac9f932922
md"""

## Simulation with variable population size

"""

# ╔═╡ 5474118e-12dc-11eb-25d2-9df31c0304cd
d

# ╔═╡ 56d2c210-12dc-11eb-3bc6-ada36466cb3a
replicate(d)

# ╔═╡ 59bf0c90-12dc-11eb-1bca-9d27f3eb8f8e
mating_extended(d, fitnesses)

# ╔═╡ bfe6dbb2-2abf-11eb-0579-31c8a0aa0bf8


# ╔═╡ 0a888190-12de-11eb-0da6-c9e35bccea74
md"""
Population size N grows to carrying capacity K and fluctuates around it when z ≈ θ. 
"""

# ╔═╡ 1f8abf20-12d6-11eb-381c-5f57e823b01c
md"""

## Mutation

"""

# ╔═╡ 2f6f3290-12d6-11eb-336a-8f3a8d7f1fb9
function mutate(a::Agent{T}, μ) where T
    newloci = Vector{T}(undef, length(a))
    for i in 1:length(a)
    if rand() > μ
            newloci[i] = a[i]
        else
            a[i] == 0.0 ? x = α : x = 0.0
            newloci[i] = x
        end
    end
    Agent(newloci)
end

# ╔═╡ 9f1bba12-12e4-11eb-2969-e5430d1d5acf
function mutation_rate(a::Agent,b::Agent)
	sum(abs.(a.loci .- b.loci)/α)/length(a)
end
	

# ╔═╡ 894e3f90-1655-11eb-3a5a-ffe9bd36d3cc
function mutate(d::Deme{A},μ) where A
    newdeme = Vector{A}(undef, length(d))
    for i=1:length(d)
		newdeme[i] = mutate(d.agents[i],μ)
	end
	Deme(newdeme, d.K, d.θ)
end 

# ╔═╡ 9c0e46b0-125a-11eb-0a79-1f336e865689
function evolving_deme_popvar_old(d::Deme, ngen, rm, Vs, μ; fun=heterozygosities, fit=malthusian_fitness, trait_mean = trait_mean, allelefreq = allelefreq)
	stats = [fun(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreq(d)]
	for n=1:ngen
		d = replicate(d)
		d = random_mating(d)
		d = mutate(d,μ)
		push!(stats, fun(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreq(d))
	end
	(stats=stats, pop = pop, tm = tm, af = af, deme=d, ngen=ngen)
end

# ╔═╡ 506821c0-207b-11eb-2bab-ddb764410ac0
function evolving_deme_popvar(d::Deme, ngen, rm, Vs, μ; fun=heterozygosities, fit=malthusian_fitness, trait_mean = trait_mean, allelefreq = allelefreq)
	stats = [fun(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreq(d)]
	for n=1:ngen
		d = mating_PnB(d,rm,Vs)
		d = mutate(d,μ)
		push!(stats, fun(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreq(d))
	end
	(stats=stats, pop = pop, tm = tm, af = af, deme=d, ngen=ngen)
end

# ╔═╡ d7bbd560-125a-11eb-14e4-8f0f82421b69
sim_popvar = evolving_deme_popvar(d, t, rm, Vs, μ)

# ╔═╡ e16d6190-12d8-11eb-267b-c55ef9f6b415
begin
	Hₒ_popvar = map(mean, sim_popvar.stats)
	plot(Hₒ_popvar, grid=false, color=:black, label="\$H_o(t)\$")
	plot!(0:sim_extended.ngen+1, 
		t->expected_heterozygosity(p*(1-p), t, N),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 6a4e2b40-2abf-11eb-181c-f3503bec0926
begin
	af_popvar = sim_popvar.af
	scatter(sort(af_popvar[50]), grid=false, color=:black, label="Random mating")
	xlabel!("Locus")
	ylabel!("Allele frequency after t generations")
end

# ╔═╡ 6cae5990-12dd-11eb-0f7d-f3d75a1710f0
begin
	popsize = map(mean, sim_popvar.pop)
	plot(popsize, grid=false, color=:black, label="Simulated population size")
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ a707c190-12eb-11eb-3caf-5d6c88a0ab92
begin
	traitmean_pv = map(mean, sim_popvar.tm)
	plot(traitmean_pv, grid=false, color=:black, label="Stabilizing selection")
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ f7d402b0-1655-11eb-1b68-ddfc8f11b202
d.agents

# ╔═╡ d7f6f030-1658-11eb-315d-5156b938eaba
md"""
## Linear habitat with migration
"""

# ╔═╡ f1948e30-1658-11eb-28ca-630b90938e3d
function random_walk(h::Habitat, p)
    new_h = emptycopy(h)
    for (i, deme) in enumerate(h.demes)
        for agent in deme.agents
            step = rand() < p ? rand([-1,1]) : 0 
            if step == -1 && i == 1
                step = 0
            elseif step == 1  && i == length(h)
                step = 0
            end
            push!(new_h.demes[i+step], agent)
        end
    end
    new_h
end

# ╔═╡ f2cfe500-2077-11eb-1103-6f128af0a6ad
function Gaussian_dispersal(h::Habitat,σ)
    new_h = emptycopy(h)
	dist = Normal(0,σ)
	dist_trunc = truncated(dist,-2*σ,2*σ)
	bin_1 = pdf(dist, σ)
	for (i, deme) in enumerate(h.demes)
        for agent in deme.agents
            step = -bin_1 < rand(dist_trunc) < bin_1  ?  0 : rand([-1,1])
            if step == -1 && i == 1
                step = 0
            elseif step == 1  && i == length(h)
                step = 0
            end
            push!(new_h.demes[i+step], agent)
        end
    end
    new_h
end

# ╔═╡ acdcfc50-2ea4-11eb-18e3-af0321250ebc
md"""
### Initiation of population
"""

# ╔═╡ dc3037a0-1659-11eb-1810-5382c488bc14
function linear_gradient(Dm,b)
    #KK = [i*b for i in 0:Dm-1]
	#KK = [θ for i in 0:Dm-1]
	a = -(Dm/2)*b + θ 
	KK = [(i*b)+a for i in 0:Dm-1]
    return KK
end

# ╔═╡ 88b085a0-2d22-11eb-36b6-4d270f2fd5df
function quadratic_gradient(Dm)
	KK = [(1/500)*(i-125)^2+12.5 for i in 0:Dm-1]
	return KK
end

# ╔═╡ e0ef0910-2d1f-11eb-0485-834f1a24426b
function cubic_gradient(Dm)
	x1 = 223/7812500
	x2 = -671/62500
	KK = [0.000499*(i^3)-0.0894*(i^2)+4.1*i for i in 0:Dm-1]
	return KK
end

# ╔═╡ 2cd1751e-2d20-11eb-256a-613452c71878
g_quad = quadratic_gradient(Dm)

# ╔═╡ dde40220-1659-11eb-04e5-bfbde2b81208
g_lin = linear_gradient(Dm,b)

# ╔═╡ 04277ec0-2d21-11eb-0a40-c71b3333b788
begin
plot(g_quad)
plot!(g_lin)
end

# ╔═╡ 73c3eda0-165a-11eb-3356-9b27c2dba452
function initiate_habitat(gradient,b,nd_s)
	"""The aim should be to initiate a population for nd_s demes for a linear gradient (with slope b) and with "optimal genetic variance" and "one half of the genes are adapted, meaning their clines take the form and spacing as assumed for the deterministic model under linkage equilibrium"""
	hab = Habitat([Deme(randagent(DiscreteNonParametric([0,α],[p,1-p]), 0, 0),K,i) for i in gradient])
	
	for j in 1:nd_s
		#x = round(Int64,Dm/2)+round(Int64,nd_s/2)-j
		#L_adapted = round(Int64,(b*x)/(2*α)) 
		
		
		for i in 1:(N)
			push!(hab.demes[round(Int64,Dm/2)+round(Int64,nd_s/2)-j].agents,randagent(DiscreteNonParametric([0,α], [p,1-p]), round(Int64,L/2)-nd_s-78)) 
			for k in 1:(round(Int64,L/2)-j-nd_s+24)
				push!(hab.demes[round(Int64,Dm/2)+12-j].agents[i], α)
			end
			while length(hab.demes[round(Int64,Dm/2)+12-j].agents[i]) < L
				push!(hab.demes[round(Int64,Dm/2)+12-j].agents[i], α-α)
			
			end
		end
	end
	hab
end

# ╔═╡ f47f5000-2d2b-11eb-2b47-055dcccaf248
hab = initiate_habitat(g_lin, b, 25)

# ╔═╡ d34571b0-2ea4-11eb-1c44-dfef7bc82957
md"""
### Simulation for linear habitat
"""

# ╔═╡ 63ef89e0-16d5-11eb-3b67-914ea2874f00
function evolving_habitat_old(h::Habitat{D}, ngen, rm, Vs, μ, p) where D
	for n = 1:ngen
		#h = random_walk(h,p)
		h = Gaussian_dispersial(h,σ)
		new_h = Vector{D}(undef, length(h))
		for (i, d) in enumerate(h.demes)
			d = replicate(d)
			d = random_mating(d)
			d = mutate(d,μ)
			new_h[i] = d
		end
		h = Habitat(new_h)
	end
	(h=h, ngen=ngen)
end

# ╔═╡ fe9b85c0-2080-11eb-069e-c730f12ff71c
function evolving_habitat(h::Habitat{D}, ngen, rm, Vs, μ, p) where D
	for n = 1:ngen
		#h = random_walk(h,p)
		h = Gaussian_dispersal(h,σ)
		new_h = Vector{D}(undef, length(h))
		for (i, d) in enumerate(h.demes)
			d = mating_PnB(d,rm,Vs)
			d = mutate(d,μ)
			new_h[i] = d
		end
		h = Habitat(new_h)
	end
	(h=h, ngen=ngen)
end

# ╔═╡ dc134e70-16d5-11eb-2c19-8dc7a5d152b7
sim_hab = evolving_habitat(hab,500,1.06,0.5,10^-6,0.50)

# ╔═╡ e30d8380-2ea4-11eb-3971-91192814178e
md"""
### Some (useful) plots
"""

# ╔═╡ 26f50f4e-17c6-11eb-336f-fdb8ab47ab73
pop_sizes = [length(deme) for deme  in sim_hab[1].demes]

# ╔═╡ 02aac282-229c-11eb-2123-e58f1f0da0e6
s = α^2/(2*Vs)

# ╔═╡ a5033c20-229b-11eb-0ee4-555c034050f5
margin = (sqrt(2)*b*σ)/((2*rm*sqrt(Vs))-b*σ) .>= 0.15.*pop_sizes.*σ*sqrt(s)

# ╔═╡ 6e5a0990-17c6-11eb-03cc-1bc05b8c5bf8
begin
	pop_sizes_p = map(mean, pop_sizes)
	p1 = plot(pop_sizes_p, grid=false, color=:black, label=false)
	hline!([K], label = "K")
	hline!([K*(1-(σ*b)*(1/(2*sqrt(Vs)*rm)))], label = "Expected pop size")
	vline!([Dm/2], label = "Starting deme")
	plot!([margin]*10, label = "Deterministic range margin")
	xlabel!("Space")
	ylabel!("Population size N")
end

# ╔═╡ c498e510-17c6-11eb-3299-a5c233e28518
trait_means = [trait_mean(deme) for deme in sim_hab[1].demes]

# ╔═╡ 8348b390-17c8-11eb-1829-f3fad57d7a02
function f_trait_agents(sim_hab)
	trait_agents = []
	cordst = []
	for (i, deme) in enumerate(sim_hab[1].demes)
		for agent in deme.agents
			t = sum(agent)
			p = (i,t)
			push!(cordst,i)
			push!(trait_agents,t)
		end
	end
	trait_agents, cordst
end

# ╔═╡ 923ba4e0-2d2e-11eb-17e3-f34323e735e8
trait_agents, cordst = f_trait_agents(sim_hab)

# ╔═╡ 5cf39490-17c7-11eb-125a-1549b47468ba
begin
	trait_means_p = map(mean, trait_means)
	p2 = plot(trait_means_p, grid=false, color=:black, label="Z_mean deme")
	plot!(g_lin, grid=false, color=:blue, label="Z optimum", linestyle=:dash)
	plot!(cordst,trait_agents, label="Z agents")
	xlabel!("Space")
	ylabel!("Trait Z")
end

# ╔═╡ fd9e9dc0-17f1-11eb-0537-63d90fb3c203
function f_het_demes(sim_hab)
	het_demes = []
	cordsh = []
	for (i, deme) in enumerate(sim_hab[1].demes)
	if length(deme) != 0
		het = α^2*sum(heterozygosities(deme))
		push!(cordsh,i)
		push!(het_demes,het)
	else
		push!(cordsh,i)
		push!(het_demes,0)
		end
	end
	het_demes, cordsh
end

# ╔═╡ 52f3b160-2d2e-11eb-3b4e-f9c40e5b8c4c
het_demes, cordsh = f_het_demes(sim_hab)


# ╔═╡ 5497f602-17fa-11eb-1fae-fbd03f58d6d6
begin
	#het_means_p = map(mean, het_demes)
	p3 = plot(cordsh, het_demes, grid=false, color=:black, label="Z_mean deme")
	hline!([b*σ*sqrt(Vs)], label = "E(V_G)")
	xlabel!("Space")
	ylabel!("\$V_G\$")
end

# ╔═╡ 8fc22ec0-27af-11eb-2977-fd6fccaef17b
begin
	
	p4 = plot(legend = false, label = "allele frequency per locus")
	xlabel!("Locus")
	ylabel!("\$p_i\$")
	for i in 1:5:length(sim_hab[1].demes)
		if length(sim_hab[1].demes[i]) != 0
		plot!(sort(allelefreq(sim_hab[1].demes[i])),linestyle=:solid, color=:black)
		end
	end
	p4

	
end

# ╔═╡ 2e160f90-2dcf-11eb-00b2-f585b00a7704
begin
hab_loci = []
for i in 1:length(sim_hab[1].demes)
	#if length(sim_hab[1].demes[i]) != 0
		for j in (sim_hab[1].demes[i].agents)
				push!(hab_loci,j.loci)
		end
	#end
end
end

# ╔═╡ 44db68d0-2de6-11eb-04b2-a19791c361da
begin
hab_het = []
for i in 1:length(sim_hab[1].demes)
	deme = (sim_hab[1].demes[i])
	if length(deme) != 0
		het = allelefreq(deme)
		push!(hab_het,het)
	else
		push!(hab_het,zeros(L))
	end
end
end

# ╔═╡ ae2ba2d0-2de8-11eb-080a-e5859a8982d7
begin
clines = []
for j in 1:Dm
	A = []
	for i in 1:L
		push!(A,(hab_het)[i][j])
	end
	push!(clines,A)
end
end

# ╔═╡ dfd18e40-2de7-11eb-081b-017c358bc0c9
sort(hab_het[90])

# ╔═╡ 538fcb40-2de7-11eb-0eb2-15bee3944ae4
begin
	p5 = plot(legend = false, label = "Clines")
	xlabel!("Space")
	ylabel!("\$p_i\$")
	for i in 1:5:length(hab_het)
		plot!((clines[i]),linestyle=:solid, color=:black)
		
	end
	p5
	
end

# ╔═╡ 8fdab3c0-17e0-11eb-21ad-0593b5d1837f
plot(p1,p2,p3,p4,p5, legend = false)

# ╔═╡ 77e8bb20-2ea3-11eb-27c6-a7cf3e978ba4
md"""
#### Simulations cfr. PnB figure 2
"""

# ╔═╡ 850a6eb2-2deb-11eb-34b7-b10e49189329
function expansion_sim(n,gen)
	"""n: number of simulated data points, gen: number of generations for each point"""
	#need to collect B and Nσ*sqrt(s)
	pop_sizes = []
	pop_diffs = []
	gen_vars = []
	B_sims = []
	SD_sims = []
	for i in 1:n
		"""Draw random parameters from range used in paper"""
		b_s = rand(Uniform(0.01,1.99))
		σ_s = rand(Uniform(0.5,4.8))
		Vs_s = rand(Uniform(0.006,8.4))
		μ_s = rand(Uniform(10^-8,10^-5))
		rm_s = rand(Uniform(0.27,2))
		K_s = rand(Uniform(4,185))
		
		"""Initiate habitat"""
		
		g = linear_gradient(Dm,b_s)
		hab = initiate_habitat(g, b_s, 25)
		sim_hab = evolving_habitat(hab,1,rm_s,Vs_s,μ_s,σ_s)
			
		"""Simulate 500 generations"""
		#sim_hab = sim_hab
		
		sim_hab_e = evolving_habitat(sim_hab[1],gen,rm_s,Vs_s,μ_s,σ_s)
		
		pop_size = sum([length(deme) for deme  in sim_hab[1].demes])
		pop_size_e = sum([length(deme) for deme  in sim_hab_e[1].demes])
		
		pop_diff = pop_size_e/pop_size
		
		gen_var = sum(f_het_demes(sim_hab)[1])
		
		push!(pop_sizes,pop_size/Dm)
		push!(gen_vars,gen_var/Dm)
		push!(pop_diffs, pop_diff)
		
		#N_s = pop_size /Dm
		N_central = K_s
		#N_e = (K_s*(1-(σ_s*b_s)))/((2*sqrt(Vs_s)*rm_s))
		s_s = (α^2)/(2*Vs_s)
		SD_s = K_s*σ_s*sqrt(s_s) 
		Vg_s = gen_var/Dm
		B_s = (b_s*σ_s)/((rm_s-(Vg_s/2*Vs_s))*sqrt(2*Vs_s))
		#B_s = (b_s*σ_s)/(((rm_s-((b_s*σ_s*sqrt(Vs_s))/2*Vs_s)))*sqrt(2*Vs_s))
		push!(B_sims, B_s)
		push!(SD_sims, SD_s)
	end
	pop_sizes, gen_vars, B_sims, SD_sims, pop_diffs
end
	

# ╔═╡ f706af60-2df5-11eb-3a8f-d53b0fbd0456
"""This line simulates data for fig. 2"""
# PS, GS, BS, SS, PD = expansion_sim(500,50)

# ╔═╡ e1d305c0-2df6-11eb-1e2c-67d412f00606
simulated_points = [(a,b) = (SS[i],BS[i]) for i in 1:length(SS)]

# ╔═╡ 0e6c9060-2df7-11eb-2c34-43b3bf3cc99c
# Note, Polechova & Barton have double logarithmic axes!
begin
	x = 0:0.1:60
	pBS = plot(x,0.15*x, legend = false, xlim=(0,60), ylim=(0,6), xscale=:log10, yscale=:log10)
	xlabel!("N*σ*sqrt(s)")
	ylabel!("B")
    colors = [x < 1. ? :red : :black for x in simulated_points]
    scatter(simulated_points, color=colors)  # I think this should work too, not tested
	#for (i,point) in enumerate(simulated_points)
	#		if PD[i] < 1
	#			scatter!([point], color=:red, label = "expansion")
	#		else PD[i] > 1
	#			scatter!([point], color=:black, label = "contraction")
	#		end
	#end
	pBS
end

# ╔═╡ 32aa4a60-2ea3-11eb-2daa-fffd0d23341f
md"""
##### What follows is a premature attempt at replicating fig. 4 from PnB that is best left ignored
"""

# ╔═╡ 107d41a0-2d2b-11eb-110f-efc8d807105f
begin
	total_gen_var = []
	total_pop_size = []
	for b = [0.025*i for i = 1:50]
		g = linear_gradient(Dm,b)
		hab = initiate_habitat(g,b,25)
		sim_hab = evolving_habitat(hab,100,rm,Vs,10^-6,0.50)
		pop_size = sum([length(deme) for deme  in sim_hab[1].demes])
		gen_var = sum(f_het_demes(sim_hab)[1])
		push!(total_pop_size, pop_size)
		push!(total_gen_var, gen_var)
	end
end

# ╔═╡ 325985e0-2d30-11eb-05cf-2b0bfc19adf5
total_gen_var

# ╔═╡ b86a3070-2d2c-11eb-0c88-8376e0674ba0
total_pop_size

# ╔═╡ 7ae17c50-2d30-11eb-0aeb-6703f543d5fa
Am = (1/10).*[Vg/((rm-(Vg/2*Vs))*Vs) for Vg in total_gen_var]

# ╔═╡ 0c9d2470-2d2f-11eb-24c8-853925635db4
Bm = [(b*σ)/((rm-(Vg/2*Vs))*sqrt(2*Vs)) for Vg in total_gen_var]

# ╔═╡ 5278c2c0-2d33-11eb-22fc-09e9ab9349eb
AB = [(a,b) = (Bm[i],Am[i]) for i in 1:length(Bm)]

# ╔═╡ 7c940e50-2d30-11eb-238c-bdebb3c2431f
begin
	Ap(B) = 10*B
	Ac(B) = (B.^2)/2
	Ae(B) = (B.*sqrt(2)).-2
	B = 0:0.1:2
	pAB = plot(B,Ap(B), legend = false)
	plot!(B,Ac(B))
	plot!(B,Ae(B))
	for point in AB
		scatter!([point])
	end
	pAB
end

# ╔═╡ Cell order:
# ╠═23dab282-10cc-11eb-0ade-c7f191741071
# ╠═020b1160-22ef-11eb-2974-e9d551ada921
# ╟─b93e3ee0-10c8-11eb-3db9-0f2cfebf6f7a
# ╠═6e53a382-10cc-11eb-0873-adb1d891562b
# ╠═1eb68d20-17e2-11eb-3166-d7bff64901ed
# ╠═29c22850-17e2-11eb-06ed-492978d8d4ea
# ╠═3fe1c320-17e2-11eb-0759-d508c9753c04
# ╠═4093d510-17e2-11eb-2157-c5a8d98967c9
# ╠═55fe6b60-10c9-11eb-18a6-e9963284f863
# ╠═48b0da00-10ca-11eb-2684-792d22df5a69
# ╠═6d3d5040-1658-11eb-0e06-095bb2d6b3a5
# ╠═17c43700-10cc-11eb-277b-bf13ce5f3505
# ╠═9920bd70-10ca-11eb-343a-4b068e4a4f16
# ╠═e8008b00-10ca-11eb-135b-35d961573890
# ╠═ccb28910-12e7-11eb-0869-8be87d5a4ae6
# ╠═6ff28710-12e8-11eb-2708-c90089167cd6
# ╠═070648f2-10cb-11eb-2f40-3573a65cb627
# ╟─7a674f20-3082-11eb-2956-01e2f8acb4be
# ╠═30074380-10cb-11eb-3600-f58f98f07fcd
# ╠═3be1ee80-10cb-11eb-275b-b73888bb0061
# ╠═f2877110-1223-11eb-1140-afde4ce678cd
# ╠═37ce3f12-124c-11eb-1cc0-05e92eaec860
# ╠═1f5a1afe-124e-11eb-0819-b7a13a875f45
# ╠═70c6d500-124e-11eb-3da9-e1bd16cc0a74
# ╠═74004a90-124d-11eb-1f21-39c3d9aac277
# ╠═e8641300-1223-11eb-1d39-091974fefa58
# ╠═e76d7e00-1223-11eb-05c8-d71557c58345
# ╠═f71a1890-1223-11eb-1b02-4980890c068d
# ╠═b1739730-1223-11eb-28ff-f56636bc7119
# ╠═f2213480-10cb-11eb-3243-59823606babf
# ╠═b3c37b90-12dc-11eb-1628-6d4571848222
# ╠═b5933190-1223-11eb-091e-9dd0af287609
# ╠═c095cc12-1223-11eb-3927-019755ca8d18
# ╠═839c8dfe-12e8-11eb-0a77-199fe2641df1
# ╠═19191fe0-2abb-11eb-1f81-d591e973be82
# ╟─8f9ee77e-1233-11eb-25f9-df1d20c7e554
# ╠═d07e5d80-10cb-11eb-2e72-3b1acc3996a1
# ╠═41e53d00-10cb-11eb-3153-c1a586773c80
# ╠═18cf418e-10d0-11eb-36b1-1d5426742a23
# ╠═0451e510-10d0-11eb-2d78-25fdaa8e5f5e
# ╟─93b4d780-10cb-11eb-06f1-556745d3c323
# ╠═8e06fc00-10cb-11eb-22d6-753424c41d8c
# ╠═bf49d880-10cd-11eb-28b4-4111e7c69a54
# ╠═24b4a580-10d1-11eb-1bcb-d775fb8b26b9
# ╠═a8fd1640-1259-11eb-259e-fdbd74535c88
# ╠═a13eaa3e-1259-11eb-2025-11cb34cf71f5
# ╠═cab06690-12de-11eb-07ae-232de34b05c3
# ╟─5017abf0-10d1-11eb-0f13-41564d34d322
# ╠═521567b0-1245-11eb-1c3a-e31e283d59c7
# ╠═873f5780-1249-11eb-105c-d16a6160ec4f
# ╠═66e5d820-1243-11eb-2162-755e424f3c4c
# ╠═447b6040-2079-11eb-0d3c-f7c17caae8dc
# ╠═1ec87310-1232-11eb-0bbc-895d245932fc
# ╠═bbe52202-1233-11eb-06d4-7f58dcfed50e
# ╠═c6b23d80-1247-11eb-2ef5-a7295f4aff4f
# ╟─052f6f12-12d9-11eb-2c06-5fac9f932922
# ╠═5474118e-12dc-11eb-25d2-9df31c0304cd
# ╠═56d2c210-12dc-11eb-3bc6-ada36466cb3a
# ╠═59bf0c90-12dc-11eb-1bca-9d27f3eb8f8e
# ╠═9c0e46b0-125a-11eb-0a79-1f336e865689
# ╠═506821c0-207b-11eb-2bab-ddb764410ac0
# ╠═d7bbd560-125a-11eb-14e4-8f0f82421b69
# ╠═e16d6190-12d8-11eb-267b-c55ef9f6b415
# ╠═6a4e2b40-2abf-11eb-181c-f3503bec0926
# ╠═bfe6dbb2-2abf-11eb-0579-31c8a0aa0bf8
# ╠═0a888190-12de-11eb-0da6-c9e35bccea74
# ╠═6cae5990-12dd-11eb-0f7d-f3d75a1710f0
# ╠═a707c190-12eb-11eb-3caf-5d6c88a0ab92
# ╟─1f8abf20-12d6-11eb-381c-5f57e823b01c
# ╠═2f6f3290-12d6-11eb-336a-8f3a8d7f1fb9
# ╠═9f1bba12-12e4-11eb-2969-e5430d1d5acf
# ╠═894e3f90-1655-11eb-3a5a-ffe9bd36d3cc
# ╠═f7d402b0-1655-11eb-1b68-ddfc8f11b202
# ╠═d7f6f030-1658-11eb-315d-5156b938eaba
# ╠═f1948e30-1658-11eb-28ca-630b90938e3d
# ╠═f2cfe500-2077-11eb-1103-6f128af0a6ad
# ╠═acdcfc50-2ea4-11eb-18e3-af0321250ebc
# ╠═dc3037a0-1659-11eb-1810-5382c488bc14
# ╠═88b085a0-2d22-11eb-36b6-4d270f2fd5df
# ╠═e0ef0910-2d1f-11eb-0485-834f1a24426b
# ╠═2cd1751e-2d20-11eb-256a-613452c71878
# ╠═04277ec0-2d21-11eb-0a40-c71b3333b788
# ╠═dde40220-1659-11eb-04e5-bfbde2b81208
# ╠═73c3eda0-165a-11eb-3356-9b27c2dba452
# ╠═f47f5000-2d2b-11eb-2b47-055dcccaf248
# ╟─d34571b0-2ea4-11eb-1c44-dfef7bc82957
# ╠═63ef89e0-16d5-11eb-3b67-914ea2874f00
# ╠═fe9b85c0-2080-11eb-069e-c730f12ff71c
# ╠═dc134e70-16d5-11eb-2c19-8dc7a5d152b7
# ╠═e30d8380-2ea4-11eb-3971-91192814178e
# ╠═26f50f4e-17c6-11eb-336f-fdb8ab47ab73
# ╠═6e5a0990-17c6-11eb-03cc-1bc05b8c5bf8
# ╠═02aac282-229c-11eb-2123-e58f1f0da0e6
# ╠═a5033c20-229b-11eb-0ee4-555c034050f5
# ╠═c498e510-17c6-11eb-3299-a5c233e28518
# ╠═8348b390-17c8-11eb-1829-f3fad57d7a02
# ╠═923ba4e0-2d2e-11eb-17e3-f34323e735e8
# ╠═5cf39490-17c7-11eb-125a-1549b47468ba
# ╠═fd9e9dc0-17f1-11eb-0537-63d90fb3c203
# ╠═52f3b160-2d2e-11eb-3b4e-f9c40e5b8c4c
# ╠═5497f602-17fa-11eb-1fae-fbd03f58d6d6
# ╠═8fc22ec0-27af-11eb-2977-fd6fccaef17b
# ╠═2e160f90-2dcf-11eb-00b2-f585b00a7704
# ╠═44db68d0-2de6-11eb-04b2-a19791c361da
# ╠═ae2ba2d0-2de8-11eb-080a-e5859a8982d7
# ╟─dfd18e40-2de7-11eb-081b-017c358bc0c9
# ╠═538fcb40-2de7-11eb-0eb2-15bee3944ae4
# ╠═8fdab3c0-17e0-11eb-21ad-0593b5d1837f
# ╠═77e8bb20-2ea3-11eb-27c6-a7cf3e978ba4
# ╠═850a6eb2-2deb-11eb-34b7-b10e49189329
# ╠═f706af60-2df5-11eb-3a8f-d53b0fbd0456
# ╠═e1d305c0-2df6-11eb-1e2c-67d412f00606
# ╠═0e6c9060-2df7-11eb-2c34-43b3bf3cc99c
# ╠═32aa4a60-2ea3-11eb-2daa-fffd0d23341f
# ╠═7c940e50-2d30-11eb-238c-bdebb3c2431f
# ╠═107d41a0-2d2b-11eb-110f-efc8d807105f
# ╠═325985e0-2d30-11eb-05cf-2b0bfc19adf5
# ╠═b86a3070-2d2c-11eb-0c88-8376e0674ba0
# ╠═7ae17c50-2d30-11eb-0aeb-6703f543d5fa
# ╠═0c9d2470-2d2f-11eb-24c8-853925635db4
# ╠═5278c2c0-2d33-11eb-22fc-09e9ab9349eb
