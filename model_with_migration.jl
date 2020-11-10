### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

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
md"""

Global parameters:

-N $(@bind N Slider(1:100, default=100, show_value=true))

-L $(@bind L Slider(1:1000, default=100, show_value=true))

-p $(@bind p Slider(0.05:0.05:0.95, default=0.5, show_value=true)) 

-t $(@bind t Slider(1:10000, default=5000, show_value=true)) 

-Vs $(@bind Vs Slider(1:0.5:10, default=2, show_value=true))

-rm $(@bind rm Slider(1.00:0.005:1.030, default=1.025, show_value=true))

-α $(@bind α Slider(0.000:0.001:5.000, default=0.25, show_value=true))

-θ $(@bind θ Slider(0:0.5:500, default=12.5, show_value=true))

-K $(@bind K Slider(1:1:500, default=125, show_value=true))

-μ $(@bind μ Slider(0.0:0.000001:1.0, default=0.000001, show_value=true))

"""

# ╔═╡ 59b194b0-17e2-11eb-1a9f-6d0fdc23f69f
σ = sqrt(1/2)

# ╔═╡ 7bad6120-17e2-11eb-2540-1fec12fbe811
h = 1

# ╔═╡ d21a5c10-17e3-11eb-0609-55c20e102990
Dm = 250

# ╔═╡ 440da79e-17e4-11eb-06b8-7f3022339dfe
b = 0.1

# ╔═╡ 1eb68d20-17e2-11eb-3166-d7bff64901ed
begin
	num_demes(σ,Vs,α) = round(Int64,10*4*sqrt(σ^2*Vs)/α)
	num_loci(nd,b,α) = round(Int64,nd*b/α)
	num_ind(K,b,σ,h,Vs,rm) = round(Int64,K*((1-b*σ)/(2*h^2*sqrt(Vs)*rm)))
end

# ╔═╡ 29c22850-17e2-11eb-06ed-492978d8d4ea
nd = num_demes(σ,Vs,α)

# ╔═╡ 3fe1c320-17e2-11eb-0759-d508c9753c04
nl = num_loci(nd,b,α)

# ╔═╡ 4093d510-17e2-11eb-2157-c5a8d98967c9
ni = num_ind(K,b,σ,h,Vs,rm)

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

# ╔═╡ 30074380-10cb-11eb-3600-f58f98f07fcd
function mate(a::Agent{T}, b::Agent{T}) where T
    newloci = Vector{T}(undef, length(a))
    for i in 1:length(a)
		newloci[i] = rand() < 0.5 ? a[i] : b[i] 
    end
	Agent(newloci)
end 

# ╔═╡ 3be1ee80-10cb-11eb-275b-b73888bb0061
function random_mating(d::Deme{A}) where {A}
    newdeme = Vector{A}(undef, length(d))
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
        	if i.loci[j] == α
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
function neutral_evolving_deme(d::Deme, ngen; fun=heterozygosities, trait_mean= trait_mean)
	stats = [fun(d)]
	tm = [trait_mean(d)]
	for n=1:ngen
		d = random_mating(d)
		push!(stats, fun(d))
		push!(tm, trait_mean(d))
	end
	(stats=stats, tm = tm, deme=d, ngen=ngen)
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
function mating_PnB(d::Deme{A}, rm, Vs) where A
	newdeme = Vector{A}(undef,0)
	fitnesses = exp.(malthusian_fitness(d,rm,Vs))
	for i=1:length(d)
		B1 = d.agents[i]
		noff = number_of_offspring(d,B1,rm,Vs)
		B2 = sample(d.agents, weights(fitnesses))
		#B2 = rand(d.agents)
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
		d = mating_extended(d, fitnesses)
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
function evolving_deme_popvar_old(d::Deme, ngen, rm, Vs, μ; fun=heterozygosities, fit=malthusian_fitness, trait_mean = trait_mean)
	stats = [fun(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	for n=1:ngen
		d = replicate(d)
		d = random_mating(d)
		d = mutate(d,μ)
		push!(stats, fun(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
	end
	(stats=stats, pop = pop, tm = tm, deme=d, ngen=ngen)
end

# ╔═╡ 506821c0-207b-11eb-2bab-ddb764410ac0
function evolving_deme_popvar(d::Deme, ngen, rm, Vs, μ; fun=heterozygosities, fit=malthusian_fitness, trait_mean = trait_mean)
	stats = [fun(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	for n=1:ngen
		d = mating_PnB(d,rm,Vs)
		d = mutate(d,μ)
		push!(stats, fun(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
	end
	(stats=stats, pop = pop, tm = tm, deme=d, ngen=ngen)
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
function Gaussian_dispersion_kernel(h::Habitat,σ)
    new_h = emptycopy(h)
	dist = Normal(0,σ)
	dist_trunc = truncated(dist,-2*σ,2*σ)
	bin_1 = pdf(dist, σ)
	for (i, deme) in enumerate(h.demes)
        for agent in deme.agents
            step = -bin_1 < rand(dist_trunc) < bin_1  ? rand([-1,1]) : 0 
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

# ╔═╡ dc3037a0-1659-11eb-1810-5382c488bc14
function linear_gradient(Dm,b)
    KK = [i*b for i in 0:Dm-1]
    return KK
end

# ╔═╡ dde40220-1659-11eb-04e5-bfbde2b81208
g = linear_gradient(Dm,b)

# ╔═╡ 73c3eda0-165a-11eb-3356-9b27c2dba452
begin
	hab = Habitat([Deme(randagent(DiscreteNonParametric([0,α],[p,1-p]), 0, 0),K,i) for i in g])
	for i in 1:N
	push!(hab.demes[round(Int64,Dm/2)].agents,randagent(DiscreteNonParametric([0,α],[p,1-p]), L))
	end
end

# ╔═╡ 63ef89e0-16d5-11eb-3b67-914ea2874f00
function evolving_habitat_old(h::Habitat{D}, ngen, rm, Vs, μ, p) where D
	for n = 1:ngen
		h = random_walk(h,p)
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
		h = Gaussian_dispersion_kernel(h,σ)
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
sim_hab = evolving_habitat(hab,5000,1.025,2,10^-6,0.50)

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

# ╔═╡ a8fdf1a0-17e0-11eb-1f45-7ddb59d9ec2d


# ╔═╡ 8348b390-17c8-11eb-1829-f3fad57d7a02
begin
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
end

# ╔═╡ 5cf39490-17c7-11eb-125a-1549b47468ba
begin
	trait_means_p = map(mean, trait_means)
	p2 = plot(trait_means_p, grid=false, color=:black, label="Z_mean deme")
	plot!(g, grid=false, color=:blue, label="Z optimum", linestyle=:dash)
	plot!(cordst,trait_agents, label="Z agents")
	xlabel!("Space")
	ylabel!("Trait Z")
end

# ╔═╡ fd9e9dc0-17f1-11eb-0537-63d90fb3c203
begin 
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
end

# ╔═╡ 5497f602-17fa-11eb-1fae-fbd03f58d6d6
begin
	#het_means_p = map(mean, het_demes)
	p3 = plot(cordsh, het_demes, grid=false, color=:black, label="Z_mean deme")
	hline!([b*σ*sqrt(Vs)], label = "E(V_G)")
	xlabel!("Space")
	ylabel!("\$V_G\$")
end

# ╔═╡ 8fdab3c0-17e0-11eb-21ad-0593b5d1837f
plot(p1,p2,p3, legend = false)

# ╔═╡ Cell order:
# ╠═23dab282-10cc-11eb-0ade-c7f191741071
# ╠═020b1160-22ef-11eb-2974-e9d551ada921
# ╟─b93e3ee0-10c8-11eb-3db9-0f2cfebf6f7a
# ╟─6e53a382-10cc-11eb-0873-adb1d891562b
# ╠═59b194b0-17e2-11eb-1a9f-6d0fdc23f69f
# ╠═7bad6120-17e2-11eb-2540-1fec12fbe811
# ╠═d21a5c10-17e3-11eb-0609-55c20e102990
# ╠═440da79e-17e4-11eb-06b8-7f3022339dfe
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
# ╠═dc3037a0-1659-11eb-1810-5382c488bc14
# ╠═dde40220-1659-11eb-04e5-bfbde2b81208
# ╠═73c3eda0-165a-11eb-3356-9b27c2dba452
# ╠═63ef89e0-16d5-11eb-3b67-914ea2874f00
# ╠═fe9b85c0-2080-11eb-069e-c730f12ff71c
# ╠═dc134e70-16d5-11eb-2c19-8dc7a5d152b7
# ╠═26f50f4e-17c6-11eb-336f-fdb8ab47ab73
# ╠═6e5a0990-17c6-11eb-03cc-1bc05b8c5bf8
# ╠═02aac282-229c-11eb-2123-e58f1f0da0e6
# ╠═a5033c20-229b-11eb-0ee4-555c034050f5
# ╠═c498e510-17c6-11eb-3299-a5c233e28518
# ╟─a8fdf1a0-17e0-11eb-1f45-7ddb59d9ec2d
# ╠═8348b390-17c8-11eb-1829-f3fad57d7a02
# ╠═5cf39490-17c7-11eb-125a-1549b47468ba
# ╠═fd9e9dc0-17f1-11eb-0537-63d90fb3c203
# ╠═5497f602-17fa-11eb-1fae-fbd03f58d6d6
# ╠═8fdab3c0-17e0-11eb-21ad-0593b5d1837f
