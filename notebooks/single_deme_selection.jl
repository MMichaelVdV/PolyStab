### A Pluto.jl notebook ###
# v0.12.17

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

# ╔═╡ b93e3ee0-10c8-11eb-3db9-0f2cfebf6f7a
md"""
# Single deme population genetics 

This notebook contains an extension of single deme population genetics with population regulation and stabilizing selection.

"""

# ╔═╡ 6e53a382-10cc-11eb-0873-adb1d891562b
md"""

Global parameters:

-N $(@bind N Slider(1:100, default=100, show_value=true))

-L $(@bind L Slider(1:1000, default=100, show_value=true))

-p $(@bind p Slider(0.05:0.05:0.95, default=0.5, show_value=true)) 

-t $(@bind t Slider(1:1000, default=500, show_value=true)) 

-Vs $(@bind Vs Slider(1:10, default=2, show_value=true))

-rm $(@bind rm Slider(1.00:0.005:1.030, default=1.025, show_value=true))

-α $(@bind α Slider(0.000:0.001:5.000, default=1.000, show_value=true))

-θ $(@bind θ Slider(0:1:500, default=50, show_value=true))

-K $(@bind K Slider(1:1:500, default=300, show_value=true))

-μ $(@bind μ Slider(0.0:0.000001:1.0, default=0.000001, show_value=true))

"""

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

# ╔═╡ 17c43700-10cc-11eb-277b-bf13ce5f3505
begin
	randagent(d::Distribution, n) = Agent(rand(d, n))
	randagent(d::Distribution, n, N) = [Agent(rand(d, n)) for i=1:N]
end

# ╔═╡ 9920bd70-10ca-11eb-343a-4b068e4a4f16
begin
	Base.length(d::Deme) = length(d.agents)
	Base.length(a::Agent) = length(a.loci)
	Base.show(io::IO, d::Deme{A,T}) where {A,T} = 
		write(io, "Deme{$A,$T}(N=$(length(d)))")
	Base.push!(d::Deme, a::Agent) = push!(d.agents, a)
	Base.push!(a::Agent{T}, x::T) where T = push!(a.loci, x)
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

# ╔═╡ 3eea0d00-124d-11eb-3b22-09d69e24887b
md"""

Am I correct that the function allelef at the moment only seems to work for α == 1 ?

"""

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

##### Here I first implemented a simulation with constant population size but mating probability proportional to fitness

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

## Mutation (tbc'd)

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

# ╔═╡ 5703ea90-12e4-11eb-21db-65b7a2302b24
a = d[1]

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
function evolving_deme_popvar(d::Deme, ngen, rm, Vs, μ; fun=heterozygosities, fit=malthusian_fitness, trait_mean = trait_mean)
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

# ╔═╡ 706bab82-12e4-11eb-270b-ffa33a42e131
b = mutate(a, 0.5)

# ╔═╡ c78c8290-12e4-11eb-3bba-7780a51e1365
mutation_rate(a,b)

# ╔═╡ f7d402b0-1655-11eb-1b68-ddfc8f11b202
d.agents

# ╔═╡ Cell order:
# ╠═23dab282-10cc-11eb-0ade-c7f191741071
# ╟─b93e3ee0-10c8-11eb-3db9-0f2cfebf6f7a
# ╠═6e53a382-10cc-11eb-0873-adb1d891562b
# ╠═55fe6b60-10c9-11eb-18a6-e9963284f863
# ╠═48b0da00-10ca-11eb-2684-792d22df5a69
# ╠═17c43700-10cc-11eb-277b-bf13ce5f3505
# ╠═9920bd70-10ca-11eb-343a-4b068e4a4f16
# ╠═e8008b00-10ca-11eb-135b-35d961573890
# ╠═ccb28910-12e7-11eb-0869-8be87d5a4ae6
# ╠═6ff28710-12e8-11eb-2708-c90089167cd6
# ╠═070648f2-10cb-11eb-2f40-3573a65cb627
# ╠═30074380-10cb-11eb-3600-f58f98f07fcd
# ╠═3be1ee80-10cb-11eb-275b-b73888bb0061
# ╟─3eea0d00-124d-11eb-3b22-09d69e24887b
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
# ╟─cab06690-12de-11eb-07ae-232de34b05c3
# ╟─5017abf0-10d1-11eb-0f13-41564d34d322
# ╠═521567b0-1245-11eb-1c3a-e31e283d59c7
# ╠═873f5780-1249-11eb-105c-d16a6160ec4f
# ╠═66e5d820-1243-11eb-2162-755e424f3c4c
# ╠═1ec87310-1232-11eb-0bbc-895d245932fc
# ╠═bbe52202-1233-11eb-06d4-7f58dcfed50e
# ╠═c6b23d80-1247-11eb-2ef5-a7295f4aff4f
# ╟─052f6f12-12d9-11eb-2c06-5fac9f932922
# ╠═5474118e-12dc-11eb-25d2-9df31c0304cd
# ╠═56d2c210-12dc-11eb-3bc6-ada36466cb3a
# ╠═59bf0c90-12dc-11eb-1bca-9d27f3eb8f8e
# ╠═9c0e46b0-125a-11eb-0a79-1f336e865689
# ╠═d7bbd560-125a-11eb-14e4-8f0f82421b69
# ╠═e16d6190-12d8-11eb-267b-c55ef9f6b415
# ╠═0a888190-12de-11eb-0da6-c9e35bccea74
# ╠═6cae5990-12dd-11eb-0f7d-f3d75a1710f0
# ╠═a707c190-12eb-11eb-3caf-5d6c88a0ab92
# ╟─1f8abf20-12d6-11eb-381c-5f57e823b01c
# ╠═2f6f3290-12d6-11eb-336a-8f3a8d7f1fb9
# ╠═5703ea90-12e4-11eb-21db-65b7a2302b24
# ╠═706bab82-12e4-11eb-270b-ffa33a42e131
# ╠═9f1bba12-12e4-11eb-2969-e5430d1d5acf
# ╠═c78c8290-12e4-11eb-3bba-7780a51e1365
# ╠═894e3f90-1655-11eb-3a5a-ffe9bd36d3cc
# ╠═f7d402b0-1655-11eb-1b68-ddfc8f11b202
