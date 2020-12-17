### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ ed04d932-342e-11eb-2d4d-f3b7a9f5615b
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI

# ╔═╡ be973290-40aa-11eb-3e67-e97df94477d0


# ╔═╡ a0c88d50-342e-11eb-2573-f13d51aabbe2
# Think of how we will represent polyploids. We will need to hold different
# chromosomes, not sure if a vector of vectors or a matrix is best. Also we
# will need to be able to scale allelic effects, using a parameter `d`. The
# latter only come in when we compute the phenotype.
# Should the mutation rate and unreduced gamete formation rate be individual
# level parameters? Probably they should...
"""
    Agent{T}
Should be an agent of arbitrary (allowed) ploidy level.
"""
@with_kw struct Agent{T,N}
    loci::Array{Array{T,1},N}
    d::T = 1. # allelic effect scaler for different ploidy levels
end

# ╔═╡ 41a71bb0-342f-11eb-2218-db85f9e52113
begin
	
Base.length(a::Agent) = length(a.loci[1]) #should be ok if all chromosomes are same length
	
ploidy(a::Agent) = length(a.loci)

# generalize to [i,j] indexing for different ploidy levels
Base.getindex(a::Agent, i) = a.loci[i]
	
Base.getindex(a::Agent, i, j) = a.loci[i][j] 

# generalize to different ploidy levels 'k'
randagent(p, α, n; d=1.) = Agent([rand(Bernoulli(p), n) * α], d)
randagent(p, α, n, N; d=1.) = [randagent(p, α, n, d=d) for i=1:N]

#Here I implemented it with 'k' as a vector of different ploidy levels to intiatiate a mixed ploidy population
randagent_p(p, α, n, k; d=1.) = Agent([(rand(Bernoulli(p), n) * α) for i=1:k], d)
randagent_p(p, α, n, k, N; d=1.) = [randagent_p(p, α, n, rand(k), d=d) for i=1:N]


	
end

# ╔═╡ b5e9f1fe-342e-11eb-1b42-57e70edfad28
# Note, I think it makes sense to have an `AbstractDeme` type, because we can
# easily imagine that when we want to implement more complex models (say with
# different fertility across ploidy levels or assortative mating or whatever)
# that we implement a new `Deme` type, but call the same functions like
# `random_mating` or whatever. So we can try to implement the main simulation
# of a habitat for the `AbstractDeme` type so that it works for any concrete 
# deme type that implements the necesary functions (such as `random_mating`)

abstract type AbstractDeme{A} 
end

# ╔═╡ 2b000160-342f-11eb-07bc-0f389915fa2a
# This is the initial thing we need, basic random mating but with gametes of
# different ploidy levels (we assume to freely fuse).
"""
    MixedPloidyDeme{A,T}
A single random-mating, mixed-ploidy level deme, most of the 'population
genetic' environment should be implemented at this level (drift, selection,
mutation). 
- `K` : Carrying capacity
- `θ` : Environmental optimum
- `rm`: Mean Malthusian fitness
- `Vs`: Variance of stabilizing selection
- `u` : Unreduced gamete formation rate
- `μ` : Mutation rate
"""
@with_kw struct MixedPloidyDeme{A,T} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 50
    θ ::T     = 12.5
    rm::T     = 1.06
    Vs::T     = 1/2
    u ::T     = 0.
    μ ::T     = 1e-6
end

# ╔═╡ 66d75f30-342f-11eb-2e72-fdbed6aa1ffd
@with_kw struct SimpleDeme{A} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 50
end

# ╔═╡ 50310830-342f-11eb-1cf0-ff654c427368
begin
	# construct a new deme from an old one, this uses the 'type as function' syntax
	(d::MixedPloidyDeme)(agents) = MixedPloidyDeme(agents, d.K, d.θ, d.rm, d.Vs, d.u, d.μ)
	(d::SimpleDeme)(agents) = SimpleDeme(agents, d.K)
end

# ╔═╡ 71e25e1e-342f-11eb-1d0a-359cb259d7f7
begin
Base.length(d::AbstractDeme) = length(d.agents)
Base.getindex(d::AbstractDeme, i) = d.agents[i]
Base.rand(rng::AbstractRNG, d::AbstractDeme, n) = rand(rng, d.agents, n)
Base.rand(d::AbstractDeme, n::Int) = rand(d.agents, n)
Base.push!(d::AbstractDeme, a::Agent) = push!(d.agents, a)
end

# ╔═╡ 55f7d480-3ce3-11eb-3fd5-df7669747ab6
begin
# NOTE: since we used similar, we have a guarantee that i does not go out of
# bounds so we can use @inbounds to gain a tiny speed nudge (but since this is
# a function that will be called a lot it could be worthwile) be careful with
# `@inbounds` though, because it can cause silent crashes...
function mate(a::Agent, b::Agent)
    newloci = similar(a.loci)
    for i in 1:length(a)
		@inbounds newloci[i] = rand() < 0.5 ? a[i] : b[i]
    end
	Agent(newloci, a.d)
end	
end

# ╔═╡ 8d96a720-342f-11eb-0b4a-d72d25dd3d5d
# won't work for mixed-ploidy
function allelefreqs(d::AbstractDeme)
    # assuming 0 is absence of allele
    f = j->mapreduce(i->d[i][j] != 0., +, 1:length(d))/length(d)
	frequencies = map(f, 1:length(d[1]))
end

# ╔═╡ 94add6a0-342f-11eb-0142-79bf9e0f5a88
# assumes haploid
function heterozygosities(d::AbstractDeme, freqs=allelefreq(d)) 
    map(p->p*(1-p), freqs)
end

# ╔═╡ d42d81f0-3b40-11eb-01ed-6fb3dda13484
#This is a deme with only haploids to check if loss of heterozygosity is still correct
d_p1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 45, [1], 50, d=1.))

# ╔═╡ f38bb420-40a6-11eb-2dd0-c94339baeb5f
#This is a deme with only diploids to check if loss of heterozygosity is still correct
d_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 45, [2], 50, d=1.))

# ╔═╡ 509f24c0-408f-11eb-3992-7d6c02ddf64a
#This is a mixed ploidy deme with both haploids and diploids
d_p = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 45, [1,2], 50, d=1.))

# ╔═╡ 890ed6e0-3b47-11eb-10bc-871e4cd5a4e6
expected_heterozygosity(H₀, t, N) = ((1.0-1.0/N)^t)*H₀

# ╔═╡ 25753ac0-3ce1-11eb-1fd2-afc42492c3da
ploidy(d::AbstractDeme) = [length(a.loci) for a in d.agents]

# ╔═╡ db86c0e0-408d-11eb-2da0-771b6e077c60
begin 
	
#This mating function assumes that different cytotypes are incompatible (i.e. when #individuals with different ploidy hybridize, they generate no offspring). This is in #accordance with the model from Levy (1975). Selfing is allowed without cost. The #dynamics might be different when you implement mating alowing for sterile hybrid #offspring that can compete for space (and affects the malthusian fitness). 


function mate_p(a::Agent, b::Agent)
	if ploidy(a) == ploidy(b) 
		newgenome = similar(a.loci)
    	for i in 1:ploidy(a) #this loops over the different chromosomes
			#rp_a = rand([1,ploidy(a)])
			#rp_b = rand([1,ploidy(b)])
			newloci = similar(a.loci[1])
				for j in 1:length(a) #this loops over the different loci for each chrosome
					@inbounds newloci[j] = rand() < 0.5 ? a[i][j] : b[i][j]
					#still need to randomize chromosome pairing
				end
			newgenome[i] = newloci
    	end
		return Agent(newgenome, a.d)
	else
		return 0
	end
end 

function random_mating_mixedp(d::AbstractDeme{A}) where A
    new_agents =  A[]
    for i=1:length(d)
		pair = mate_p(rand(d, 2)...)
		if pair != 0
			push!(new_agents, pair)
		end
	end
	MixedPloidyDeme(agents = new_agents)
end 

end

# ╔═╡ 7f638d30-342f-11eb-0ca3-1fbe88cb80ed
begin
# implemented for abstract deme type
function random_mating(d::AbstractDeme)
    newdeme = similar(d.agents)
    for i=1:length(d)
		@inbounds newdeme[i] = mate(rand(d, 2)...)
	end
	d(newdeme)
end 

function random_mating_p(d::AbstractDeme)
    newdeme = similar(d.agents)
    for i=1:length(d)
		@inbounds newdeme[i] = mate_p(rand(d, 2)...)
	end
	d(newdeme)
end

end

# ╔═╡ 95641070-3b44-11eb-322c-49dc6cc53863
function allelefreqs_p(d::AbstractDeme)
	freq = Vector{Float64}(undef,length(d.agents[1]))
	for loc in 1:length(freq)
		s = 0
    	for ag in d.agents
			for chr in ploidy(ag)
        		if ag.loci[chr][loc] != 0
            		s += 1
				end
			end
        end
		f = s/length(d)
		freq[loc] = f
	end
	freq
end	

# ╔═╡ 66262ee0-3b46-11eb-0f77-4f65686d8a8b
function heterozygosities_p(d::AbstractDeme, freqs=allelefreqs_p(d)) 
    map(p->p*(1-p), freqs)
end

# ╔═╡ db9c5d80-3ce4-11eb-33c8-01744ab5889f
ploidy(d_p)

# ╔═╡ 46f18430-4095-11eb-3e39-e3992cd8e281
#This is a random function I quickly wrote to check the population size of haploids and diploids seperately. Should be generalized.

function ploidy_freq(d::AbstractDeme)
	p1 = 0
	p2 = 0
	for agent in d.agents
		if ploidy(agent) == 1
			p1 += 1
		elseif ploidy(agent) == 2
			p2 += 1
		end
	end
	p1, p2
end		

# ╔═╡ 4e3640e0-3b46-11eb-3bac-59e50426df3c
function neutral_evolving_deme(d::AbstractDeme, ngen; heterozygosities_p = heterozygosities_p, allelefreqs_p = allelefreqs_p, trait_mean = trait_mean, pf = ploidy_freq)
	het = [heterozygosities_p(d)]
	af = [allelefreqs_p(d)]
	tm = [trait_mean(d)]
	p1 = [ploidy_freq(d)[1]]
	p2 = [ploidy_freq(d)[2]]
	
	for n=1:ngen
		d = random_mating_mixedp(d)
		push!(het, heterozygosities_p(d))
		push!(af, allelefreqs_p(d))
		push!(tm, trait_mean(d))
		push!(p1, ploidy_freq(d)[1])
		push!(p2, ploidy_freq(d)[2])
	end
	(het=het, af=af, tm=tm, deme=d, p1=p1, p2=p2, ngen=ngen)
end

# ╔═╡ 55d32ca0-40a0-11eb-1101-8b3f1beb0f2e
sim_haploid = neutral_evolving_deme(d_p1, 500)

# ╔═╡ 0e000310-40a7-11eb-05ad-0d2e3507427e
sim_diploid = neutral_evolving_deme(d_p2, 500)

# ╔═╡ 8dc157e0-3b46-11eb-3acc-6d28b136ca3c
sim = neutral_evolving_deme(d_p, 500)

# ╔═╡ 95418262-3b46-11eb-25d6-c1fe78fc90ca
begin
	p = 0.5
	t = 500
	N = 50
	Hₒ = map(mean, sim.het)
	plot(Hₒ, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH mixed")
	plot!(0:sim.ngen+1, 
		t->expected_heterozygosity(p*(1-p), t, N),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 60f31320-40a0-11eb-236c-c9933d4c36e8
begin
	Hₒ_hap = map(mean, sim_haploid.het)
	plot(Hₒ_hap, grid=false, color=:black, label="\$H_o(t)\$ ", title = "LOH haploid")
	plot!(0:sim.ngen+1, 
		t->expected_heterozygosity(p*(1-p), t, N),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
	
end

# ╔═╡ 68f64c1e-40a7-11eb-3221-57d68a17a6d9
begin
	Hₒ_dip = map(mean, sim_diploid.het)
	plot(Hₒ_dip, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH diploid")
	plot!(0:sim_diploid.ngen+1, 
		t->expected_heterozygosity(p*(1-p), t, N),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 84a53e60-3b4b-11eb-346f-a92bf779d279
begin
	af = sim.af
	scatter(af[length(af)], grid=false, color=:black, label= false)
	xlabel!("Locus")
	ylabel!("Allele frequency after t generations")
end

# ╔═╡ ea9d09b0-409f-11eb-082b-8778a053106b
#After the simulation only 1 cytotype remains in the population, also plotted this below. 
ploidy(sim.deme)

# ╔═╡ efcac440-4095-11eb-1e4e-6bf76dba2120
begin
	pf1 = sim.p1
	pf2 = sim.p2
	plot(pf1, grid=false, color=:blue, label="haploid")
	plot!(pf2, grid=false, color=:red, label="diploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ 6ca05f30-3e74-11eb-08ac-0faed2e991fa
Base.sum(a::Agent) = sum(sum(a.loci))

# ╔═╡ 3a45a2b0-3e75-11eb-2fff-2f8e200cb29d
function trait(a::Agent)
	sum(a)/a.d
end

# ╔═╡ 0f52dfb0-3e74-11eb-3f07-094774d83171
function trait_mean(d::AbstractDeme)
	z = Float64[]
	for agent in d.agents
		push!(z,trait(agent))
	end
	sum(z)/length(d)
end

# ╔═╡ 93dafc30-3e75-11eb-3f02-a7f994da33e3
begin
	traitmean = map(mean, sim.tm)
	plot(traitmean, grid=false, color=:black, label=false, title="Random mating without selection")
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ a444e532-4097-11eb-0a27-89e61edc9802
md"""
## Including population regulation and stabilizing selection
"""

# ╔═╡ 4903f02e-3e76-11eb-3d88-d1566b3a3523
function malthusian_fitness(d::AbstractDeme,a::Agent)
    N = length(d)
    z = trait(a)
    return d.rm*(1-(N/d.K))-((z-d.θ)^2)/(2*d.Vs)
end 

# ╔═╡ 4e2a4c30-3e76-11eb-0d8f-1b919ab9d04b
function malthusian_fitness(d::AbstractDeme)
    N = length(d)
    fitnesses = Float64[]
	for agent in d.agents
		z = sum(agent)
    	f = d.rm*(1-(N/d.K))-((z-d.θ)^2)/(2*d.Vs)
		push!(fitnesses, f)
	end
	fitnesses
end 

# ╔═╡ 59c7ee30-3e76-11eb-07ee-bb3b33a9f6db
number_of_offspring(d::AbstractDeme,a::Agent) = rand(Poisson(exp(malthusian_fitness(d::AbstractDeme,a::Agent))))

# ╔═╡ 7f38e42e-3e76-11eb-26b1-9d54c06c0c34
function mating_PnB(d::AbstractDeme{A}) where A
	new_agents =  A[]
	fitnesses = exp.(malthusian_fitness(d))
	for i=1:length(d)
		B1 = d.agents[i]
		noff = number_of_offspring(d,B1)
		B2 = sample(d.agents, weights(fitnesses))
		#B2 = rand(d.agents)
		#child = mate(B1,B2)
		for c in 1:noff 
				push!(new_agents, mate_p(B1,B2))
		end
	end
	MixedPloidyDeme(agents=new_agents)
end 	

# ╔═╡ 88e341b0-3e76-11eb-05d0-9d9b75e16a6a
function evolving_deme_popvar(d::AbstractDeme, ngen; heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = allelefreqs_p, pf = ploidy_freq)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p1 = [ploidy_freq(d)[1]]
	p2 = [ploidy_freq(d)[2]]
	for n=1:ngen
		d = mating_PnB(d)
		#d = mutate(d) #is this defined on level of deme or agent ?
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p1, ploidy_freq(d)[1])
		push!(p2, ploidy_freq(d)[2])
		
	end
	(het=het, pop=pop, tm=tm, af=af, deme=d, p1=p1, p2=p2, ngen=ngen)
end

# ╔═╡ 9eee9130-3e76-11eb-21e2-4516ca929e40
sim_popvar = evolving_deme_popvar(d_p, t)

# ╔═╡ a3966e60-3e76-11eb-2e07-3f8537c4a8b4
begin
	Hₒ_popvar = map(mean, sim_popvar.het)
	plot(Hₒ_popvar, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH mixed ploidy with selection")
	plot!(0:sim_popvar.ngen+1, 
		t->expected_heterozygosity(p*(1-p), t, N),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ a87f4c30-3e76-11eb-3ad3-33bf4d5d9089
begin
	af_popvar = sim_popvar.af
	scatter(sort(af_popvar[50]), grid=false, color=:black, label=false)
	xlabel!("Locus")
	ylabel!("Allele frequency after t generations")
end

# ╔═╡ aed0bab0-3e76-11eb-338f-efff9f2811e6
begin
	popsize = map(mean, sim_popvar.pop)
	plot(popsize, grid=false, color=:black, label=false, title="Simulated population size")
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ b01bb140-3e76-11eb-2716-958dbfc34606
begin
	traitmean_pv = map(mean, sim_popvar.tm)
	plot(traitmean_pv, grid=false, color=:black, label=false, title="Effect of stabilizing selection")
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ d8a5b240-4093-11eb-1b6e-f1129b1ac5b5
begin
	pf1_var = sim_popvar.p1
	pf2_var = sim_popvar.p2
	plot(pf1_var, grid=false, color=:blue, label="haploid")
	plot!(pf2_var, grid=false, color=:red, label="diploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ 3147a3b0-4097-11eb-34be-bf35f43bd3e1
md""" ## Rubbish below"""

# ╔═╡ b42552a0-3e76-11eb-3ee2-a5bf1aabcff6
function mutate(d::AbstractDeme, a::Agent) 
    newloci = similar(a.loci)
	for i in 1:ploidy(a)
    	for j in 1:length(a)
    		if rand() > d.μ
            	newloci[i][j] = a.loci[i][j]
        	else
            	a.loci[i][j] == 0.0 ? x = α : x = 0.0
            	newloci[i][j] = x
        	end
		end
    end
    Agent(newloci)
end

# ╔═╡ 619f8c50-408b-11eb-184a-95d63ec67247
AAA = d_p.agents[1]

# ╔═╡ bdd862b0-3e76-11eb-09fa-6fa07f6d6573
function mutate(d::AbstractDeme{A}) where A
    newdeme = Vector{A}(undef, length(d))
    for i=1:length(d)
		newdeme[i] = mutate(d.agents[i])
	end
	Deme(newdeme, d.K, d.θ)
end 

# ╔═╡ 9b9bac30-342f-11eb-120b-9d2e7354f241
"""
    Habitat{D}
A 1-dimensional habitat, i.e. an array of connected demes. This implements
the migration aspects of the population genetic environment.
"""
@with_kw struct Habitat{D,T}
    demes::Vector{D}
    σ ::T = 1/2
    b ::T = 0.1
    #Dm::T 
end

# ╔═╡ a36e3110-4093-11eb-026a-d3ddc3217805
emptycopy(h::Habitat) = Habitat(emptycopy.(h.demes))

# ╔═╡ ce14f100-4012-11eb-0daf-19617768225e
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

# ╔═╡ d124a200-4012-11eb-3644-2be0a62eab74
function linear_gradient(Dm,b)
    #KK = [i*b for i in 0:Dm-1]
	#KK = [θ for i in 0:Dm-1]
	a = -(Dm/2)*b + θ 
	KK = [(i*b)+a for i in 0:Dm-1]
    return KK
end

# ╔═╡ bcd073c0-4093-11eb-0484-b73fb9772521
linear_gradient(100,0.1)

# ╔═╡ d98408f0-4012-11eb-3515-bfd13ab42f56
g_lin = linear_gradient(d_p,b)

# ╔═╡ f9474a80-4012-11eb-1b5c-afa826771136
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

# ╔═╡ ff24c180-4012-11eb-3edb-09d780162b72
hab = initiate_habitat(g_lin, b, 25)

# ╔═╡ 00f20770-4013-11eb-172e-7d4cfede7a85
function evolving_habitat(h::Habitat{D}, ngen) where D
	for n = 1:ngen
		h = random_walk(h)
		#h = Gaussian_dispersal(h,σ)
		new_h = Vector{D}(undef, length(h))
		for (i, d) in enumerate(h.demes)
			d = mating_PnB(d)
			d = mutate(d)
			new_h[i] = d
		end
		h = Habitat(new_h)
	end
	(h=h, ngen=ngen)
end

# ╔═╡ 0832b7f0-4013-11eb-0d7e-d3f77e87aa9b
sim_hab = evolving_habitat(hab,1,1.06,0.5,10^-6,0.50)

# ╔═╡ 375696f0-4013-11eb-1779-677909ec6e11
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

# ╔═╡ 2e97f900-4013-11eb-316c-650ea3a40752
margin = (sqrt(2)*b*σ)/((2*rm*sqrt(Vs))-b*σ) .>= 0.15.*pop_sizes.*σ*sqrt(s)

# ╔═╡ 0f7d5380-4013-11eb-2a3b-614df9653550
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

# ╔═╡ 472b52f0-4013-11eb-1e90-4775755333ef
trait_means = [trait_mean(deme) for deme in sim_hab[1].demes]

# ╔═╡ 3ddde280-4013-11eb-39a3-130d6ec349ff
trait_agents, cordst = f_trait_agents(sim_hab)

# ╔═╡ 51427c00-4013-11eb-2be7-511a95abdf44
begin
	trait_means_p = map(mean, trait_means)
	p2 = plot(trait_means_p, grid=false, color=:black, label="Z_mean deme")
	plot!(g_lin, grid=false, color=:blue, label="Z optimum", linestyle=:dash)
	plot!(cordst,trait_agents, label="Z agents")
	xlabel!("Space")
	ylabel!("Trait Z")
end

# ╔═╡ Cell order:
# ╠═ed04d932-342e-11eb-2d4d-f3b7a9f5615b
# ╟─be973290-40aa-11eb-3e67-e97df94477d0
# ╠═a0c88d50-342e-11eb-2573-f13d51aabbe2
# ╠═41a71bb0-342f-11eb-2218-db85f9e52113
# ╠═b5e9f1fe-342e-11eb-1b42-57e70edfad28
# ╠═2b000160-342f-11eb-07bc-0f389915fa2a
# ╠═50310830-342f-11eb-1cf0-ff654c427368
# ╠═66d75f30-342f-11eb-2e72-fdbed6aa1ffd
# ╠═71e25e1e-342f-11eb-1d0a-359cb259d7f7
# ╠═55f7d480-3ce3-11eb-3fd5-df7669747ab6
# ╠═7f638d30-342f-11eb-0ca3-1fbe88cb80ed
# ╠═db86c0e0-408d-11eb-2da0-771b6e077c60
# ╠═8d96a720-342f-11eb-0b4a-d72d25dd3d5d
# ╠═95641070-3b44-11eb-322c-49dc6cc53863
# ╠═94add6a0-342f-11eb-0142-79bf9e0f5a88
# ╠═66262ee0-3b46-11eb-0f77-4f65686d8a8b
# ╠═d42d81f0-3b40-11eb-01ed-6fb3dda13484
# ╠═f38bb420-40a6-11eb-2dd0-c94339baeb5f
# ╠═509f24c0-408f-11eb-3992-7d6c02ddf64a
# ╠═4e3640e0-3b46-11eb-3bac-59e50426df3c
# ╠═55d32ca0-40a0-11eb-1101-8b3f1beb0f2e
# ╠═890ed6e0-3b47-11eb-10bc-871e4cd5a4e6
# ╠═60f31320-40a0-11eb-236c-c9933d4c36e8
# ╠═0e000310-40a7-11eb-05ad-0d2e3507427e
# ╠═68f64c1e-40a7-11eb-3221-57d68a17a6d9
# ╠═8dc157e0-3b46-11eb-3acc-6d28b136ca3c
# ╠═95418262-3b46-11eb-25d6-c1fe78fc90ca
# ╠═84a53e60-3b4b-11eb-346f-a92bf779d279
# ╠═25753ac0-3ce1-11eb-1fd2-afc42492c3da
# ╠═db9c5d80-3ce4-11eb-33c8-01744ab5889f
# ╠═ea9d09b0-409f-11eb-082b-8778a053106b
# ╠═46f18430-4095-11eb-3e39-e3992cd8e281
# ╠═efcac440-4095-11eb-1e4e-6bf76dba2120
# ╠═6ca05f30-3e74-11eb-08ac-0faed2e991fa
# ╠═3a45a2b0-3e75-11eb-2fff-2f8e200cb29d
# ╠═0f52dfb0-3e74-11eb-3f07-094774d83171
# ╠═93dafc30-3e75-11eb-3f02-a7f994da33e3
# ╠═a444e532-4097-11eb-0a27-89e61edc9802
# ╠═4903f02e-3e76-11eb-3d88-d1566b3a3523
# ╠═4e2a4c30-3e76-11eb-0d8f-1b919ab9d04b
# ╠═59c7ee30-3e76-11eb-07ee-bb3b33a9f6db
# ╠═7f38e42e-3e76-11eb-26b1-9d54c06c0c34
# ╠═88e341b0-3e76-11eb-05d0-9d9b75e16a6a
# ╠═9eee9130-3e76-11eb-21e2-4516ca929e40
# ╠═a3966e60-3e76-11eb-2e07-3f8537c4a8b4
# ╠═a87f4c30-3e76-11eb-3ad3-33bf4d5d9089
# ╠═aed0bab0-3e76-11eb-338f-efff9f2811e6
# ╠═b01bb140-3e76-11eb-2716-958dbfc34606
# ╠═d8a5b240-4093-11eb-1b6e-f1129b1ac5b5
# ╠═3147a3b0-4097-11eb-34be-bf35f43bd3e1
# ╠═b42552a0-3e76-11eb-3ee2-a5bf1aabcff6
# ╠═619f8c50-408b-11eb-184a-95d63ec67247
# ╠═bdd862b0-3e76-11eb-09fa-6fa07f6d6573
# ╠═9b9bac30-342f-11eb-120b-9d2e7354f241
# ╠═a36e3110-4093-11eb-026a-d3ddc3217805
# ╠═ce14f100-4012-11eb-0daf-19617768225e
# ╠═d124a200-4012-11eb-3644-2be0a62eab74
# ╠═bcd073c0-4093-11eb-0484-b73fb9772521
# ╠═d98408f0-4012-11eb-3515-bfd13ab42f56
# ╠═f9474a80-4012-11eb-1b5c-afa826771136
# ╠═ff24c180-4012-11eb-3edb-09d780162b72
# ╠═00f20770-4013-11eb-172e-7d4cfede7a85
# ╠═0832b7f0-4013-11eb-0d7e-d3f77e87aa9b
# ╠═0f7d5380-4013-11eb-2a3b-614df9653550
# ╠═375696f0-4013-11eb-1779-677909ec6e11
# ╠═2e97f900-4013-11eb-316c-650ea3a40752
# ╠═472b52f0-4013-11eb-1e90-4775755333ef
# ╠═3ddde280-4013-11eb-39a3-130d6ec349ff
# ╠═51427c00-4013-11eb-2be7-511a95abdf44
