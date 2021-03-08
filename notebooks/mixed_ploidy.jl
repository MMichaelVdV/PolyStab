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
    K ::Int64 = 100
    θ ::T     = 12.5
    rm::T     = 1.06
    Vs::T     = 1/2
    u ::T     = 0.01
    μ ::T     = 1e-6
end

# ╔═╡ 66d75f30-342f-11eb-2e72-fdbed6aa1ffd
@with_kw struct SimpleDeme{A} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 15
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
d_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 45, [2], 50, d=2.))

# ╔═╡ 509f24c0-408f-11eb-3992-7d6c02ddf64a
#This is a mixed ploidy deme with both diploids and tetraploids
d_p = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 0.5, 45, [2], 100, d=2.),randagent_p(0.5, 0.5, 45, [4], 0, d=4.)))

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
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
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
	p2 = 0
	p3 = 0
	p4 = 0
	for agent in d.agents
		if ploidy(agent) == 2
			p2 += 1
		elseif ploidy(agent) == 3
			p3 += 1
		elseif ploidy(agent) == 4
			p4 += 1
		end
	end
	p2, p3, p4
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

# ╔═╡ a444e532-4097-11eb-0a27-89e61edc9802
md"""
## Including population regulation and stabilizing selection
"""

# ╔═╡ 473ff8d0-4186-11eb-25a9-f71f8fdfd1d9
function unreduced_gamete(d::AbstractDeme, a::Agent)
	loci = [ Float64[] for x in 1:2*ploidy(a) ]
	if ploidy(a) == 2
		if rand() < d.u
			loci[1] = a.loci[1]
			loci[2] = a.loci[1]
			loci[3] = a.loci[2]
			loci[4] = a.loci[2]
			return Agent(loci, 4.)
		end
	end
		return a
end
		

# ╔═╡ 3637a430-4199-11eb-020d-951aeb617f11
function unreduced_gamete(d::AbstractDeme{A}) where A
	new_agents =  A[]
	for agent in d.agents
		push!(new_agents, unreduced_gamete(d,agent))
	end
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end	

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
		z = trait(agent)
    	f = d.rm*(1-(N/d.K))-((z-d.θ)^2)/(2*d.Vs)
		push!(fitnesses, f)
	end
	fitnesses
end 

# ╔═╡ 59c7ee30-3e76-11eb-07ee-bb3b33a9f6db
number_of_offspring(d::AbstractDeme,a::Agent) = rand(Poisson(exp(malthusian_fitness(d::AbstractDeme,a::Agent))))

# ╔═╡ 7d3cfed0-4184-11eb-2b59-29c1afda519d
[number_of_offspring(d_p,d_p.agents[i]) for i=1:length(d_p.agents)]

# ╔═╡ 55d3df82-729e-11eb-23b0-03eeb621d7b2


# ╔═╡ b5f6d590-717f-11eb-01dc-fb779d0a2410
md"""
### Implementation with differential offspring viability and unreduced gamete formation
"""

# ╔═╡ b4ea1900-7099-11eb-3159-453d15bc2f7a
#Viability matrix, should probably become integrated within deme struct
struct OffspringViability{T,N}
    viability::Array{Array{T,1},N}
end

# ╔═╡ b4e7f620-7099-11eb-0954-f3d8ab8d5103
#Unreduced gamete probability matrix, should probably become integrated within deme struct
struct UnreducedGamete{T,N} 
    prob::Array{Array{T,1},N}
end

# ╔═╡ a4b40320-7099-11eb-27cb-1f27ef2645ee
#Allocates some random probabalities for unreduced gamete formation for 2n, 3n, 4n individuals respectively, for example [0.9,0.1,0.,0.] -> 2n ind has probability of 0.9 to produce 1n gametes, 0.1 for 2n (unreduced) gamete,...this doesn't take into account aneuploid gametes yet.
UG = UnreducedGamete([[0.95,0.05,0.,0.],[0.,0.,0.,0.],[0.0,1.,0.0,0.]])

# ╔═╡ 940db692-70aa-11eb-1276-9d757f1ed1ba
#On the rows/columns (symmetric matrix) are the ploidy levels of the gametes to be combined (1n to 4n).The upper bound is now at 4n, all ploidy levels above are not viable (ie. combination of 2n and 3n gamete has 0 viability).
OV = OffspringViability([[1.,0.,0.1,0.],[0.,1.,0.,0.],[0.1,0.,0.,0.],[0.,0.,0.,0.]])

# ╔═╡ fdd49290-7ba1-11eb-3316-45fa9dc62000
ploidy.(d_p.agents)

# ╔═╡ 2260a980-70af-11eb-2688-4d1b7d6f2f8b
function viability(a::Agent, b::Agent, OV::OffspringViability)
	return OV.viability[ploidy(a)][ploidy(b)]
end

# ╔═╡ 23b9e2ce-70ad-11eb-1e5d-1f8ad0c849bb
function unreduced_gamete(a::Agent, UG::UnreducedGamete)
	#this samples the ploidy level (1 to 4, potentially) of gametes 
	num = sample([1.,2.,3.,4.],weights(UG.prob[ploidy(a)-1]))
	loci = [ Float64[] for x in 1:num ]
	b = deepcopy(a.loci)
	shuffle!(b)
	c = 1
	while c <= num
		i = b[1]
		popfirst!(b)
		loci[c] = i
		c += 1
	end
	return Agent(loci, num)
end	

# ╔═╡ 14e8d250-7bbd-11eb-0da7-afeb3fc81d83
function recombine_poly(a::Agent)
	num = ploidy(a)
	loci = [ Float64[] for x in 1:num ]
	al = a.loci
	newlocus = similar(al[1])
		
	for l in 1:num
		for j in 1:length(al[1]) #this loops over the different loci for each chrosome
			i = rand([x for x in 1:num])
			@inbounds newlocus[j] = al[i][j]
		end
		loci[l] = newlocus
	end
	Agent(loci=loci)
end	

# ╔═╡ 321bbf50-7bad-11eb-0b6d-6f83a6916170
function recombine_poly(d::AbstractDeme{A}) where A
	new_agents =  A[]
	for agent in d.agents
		num = ploidy(agent)
		loci = [ Float64[] for x in 1:num ]
		a = agent.loci
		newlocus = similar(a[1])
		
		for l in 1:num
			for j in 1:length(a[1]) #this loops over the different loci for each chrosome
				i = rand([x for x in 1:num])
				@inbounds newlocus[j] = a[i][j]
			end
			loci[l] = newlocus
		end
		
		push!(new_agents, Agent(loci=loci))
	end
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end		

# ╔═╡ 737a7620-7bb3-11eb-14c2-83beec620ea1
recombine_poly(d_p)

# ╔═╡ ad2350a0-7bb2-11eb-2517-1f70162a1a54
rand([x for x in 1:3])

# ╔═╡ eabf2cb0-7bb0-11eb-3e49-db4d6b9a21e2
a = [[1,2,3],[4,5,6],[7,8,9]]

# ╔═╡ 44dbde9e-7bb1-11eb-0d1b-4b44d5879262
a[2][1]

# ╔═╡ 73b40ea0-7bb1-11eb-264e-f7120f56fba4
shuffle([1,2,3])

# ╔═╡ 1d6cff90-715a-11eb-2648-c9e1455933bd
function unreduced_gamete(d::AbstractDeme{A}, UG::UnreducedGamete) where A
	new_agents =  A[]
	for agent in d.agents
		push!(new_agents, unreduced_gamete(agent,UG))
	end
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end	

# ╔═╡ 4e3640e0-3b46-11eb-3bac-59e50426df3c
function neutral_evolving_deme(d::AbstractDeme, ngen; heterozygosities_p = heterozygosities_p, allelefreqs_p = allelefreqs_p, trait_mean = trait_mean, pf = ploidy_freq)
	het = [heterozygosities_p(d)]
	af = [allelefreqs_p(d)]
	tm = [trait_mean(d)]
	p2 = [ploidy_freq(d)[1]]
	p3 = [ploidy_freq(d)[2]]
	p4 = [ploidy_freq(d)[3]]
	
	for n=1:ngen
		d = random_mating_mixedp(d)
		d = unreduced_gamete(d)
		push!(het, heterozygosities_p(d))
		push!(af, allelefreqs_p(d))
		push!(tm, trait_mean(d))
		push!(p2, ploidy_freq(d)[1])
		push!(p3, ploidy_freq(d)[2])
		push!(p4, ploidy_freq(d)[3])
	end
	(het=het, af=af, tm=tm, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen)
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
	plot!(0:sim_haploid.ngen+1, 
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
	pf1 = sim.p2
	pf2 = sim.p4
	plot(pf1, grid=false, color=:blue, label="diploid")
	plot!(pf2, grid=false, color=:red, label="tetraploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ 93dafc30-3e75-11eb-3f02-a7f994da33e3
begin
	traitmean = map(mean, sim.tm)
	plot(traitmean, grid=false, color=:black, label=false, title="Random mating without selection")
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 01be8d80-715b-11eb-3b9d-17a328be8383
begin 
	
#This mating function assumes that different cytotypes can be compatible with a decrease in viability (cfr. OffspringViability matrix) (i.e. when #individuals with different ploidy hybridize, they generate have a probability p to generate no viable offspring). Selfing is allowed without cost. This influences the dynamics by including hybrid offspring that can compete for space (and affects the malthusian fitness/density dependence of selection). Still need to incorporate a decent recombination method. 


function mate_p(a::Agent, b::Agent, UG::UnreducedGamete, OV::OffspringViability)
	#gamete formation
	ag = recombine_poly(unreduced_gamete(a,UG))
	bg = recombine_poly(unreduced_gamete(b,UG))
	gam_a = ag.loci	
	gam_b = bg.loci
	#combine gametes and assign viability
	via = viability(ag, bg, OV)
	if rand() < via
		num = ploidy(ag) + ploidy(bg)
		loci = [ Float64[] for x in 1:num ]
	
		c = 1
		while c <= ploidy(ag)
			loci[c] = gam_a[c]
			c+=1
		end
		while c <= num
			loci[c] = gam_b[c-ploidy(ag)]
			c+=1
		end
		return Agent(loci, 1. *num)
		end
	return 0
end
			
		
		
#	if ploidy(a) == ploidy(b) 
#		newgenome = similar(a.loci)
#    	for i in 1:ploidy(a) #this loops over the different chromosomes
#			#rp_a = rand([1,ploidy(a)])
#			#rp_b = rand([1,ploidy(b)])
#			newloci = similar(a.loci[1])
#				for j in 1:length(a) #this loops over the different loci for each #chrosome
#					@inbounds newloci[j] = rand() < 0.5 ? a[i][j] : b[i][j]
#					#still need to randomize chromosome pairing
#				end
#			newgenome[i] = newloci
#    	end
#		return Agent(newgenome, a.d)
#	else
#		return 0
#	end
#end 

#function random_mating_mixedp(d::AbstractDeme{A}, OV::OffspringViability) where A
#    new_agents =  A[]
#    for i=1:length(d)
#		pair = mate_p(rand(d, 2)...)
#		if pair != 0
#			push!(new_agents, pair)
#		end
#	end
#	MixedPloidyDeme(agents = new_agents)
#end 
#
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
	d(agents=newdeme,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end

end

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
		if mate_p(B1,B2) != 0
			for c in 1:noff 
				push!(new_agents, mate_p(B1,B2))
			end
		end
	end
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end 	

# ╔═╡ b0ef24ee-715a-11eb-0dd3-a9be75572bdb
function mating_PnB(d::AbstractDeme{A}, UG::UnreducedGamete, OV::OffspringViability) where A
	new_agents =  A[]
	fitnesses = exp.(malthusian_fitness(d))
	for i=1:length(d)
		B1 = d.agents[i]
		noff = number_of_offspring(d,B1)
		B2 = sample(d.agents, weights(fitnesses))
		#B2 = rand(d.agents)
		#child = mate(B1,B2)
		m = mate_p(B1,B2,UG,OV)
		if m != 0
			for c in 1:noff 
				push!(new_agents, m)
			end
		end
	end
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end 

# ╔═╡ ea898c72-40c5-11eb-233f-074fe4f50a3c
mating_PnB(d_p)

# ╔═╡ 88e341b0-3e76-11eb-05d0-9d9b75e16a6a
function evolving_deme_popvar(d::AbstractDeme, ngen; heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = allelefreqs_p, pf = ploidy_freq)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[1]]
	p3 = [ploidy_freq(d)[2]]
	p4 = [ploidy_freq(d)[3]]
	
	for n=1:ngen
		d = mating_PnB(d)
		d = unreduced_gamete(d)
		#d = mutate(d) #is this defined on level of deme or agent ?
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[1])
		push!(p3, ploidy_freq(d)[2])
		push!(p4, ploidy_freq(d)[3])
		
	end
	(het=het, pop=pop, tm=tm, af=af, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen)
end

# ╔═╡ 9eee9130-3e76-11eb-21e2-4516ca929e40
sim_popvar = evolving_deme_popvar(d_p,t)

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
	scatter(sort(af_popvar[1]), grid=false, color=:black, label=false)
	xlabel!("Locus")
	ylabel!("Allele frequency after t generations")
end

# ╔═╡ 90693580-72a0-11eb-2127-e1ca7c4c1cbc
begin
	anim = @animate for i ∈ 1:t
	    scatter(sort(af_popvar[i]), grid=false, color=:black, label=false)
	end every 1
	gif(anim, "fizzywop.gif", fps = 30)
end

# ╔═╡ 5666a9a0-729e-11eb-22ee-cfd529fe980b
af_popvar

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
	pf1_var = sim_popvar.p2
	pf2_var = sim_popvar.p4
	plot(pf1_var, grid=false, color=:blue, label="diploid")
	plot!(pf2_var, grid=false, color=:red, label="tetraploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ 313d5100-714b-11eb-0f11-db84567a7e36
function evolving_deme_ploidyvar(d::AbstractDeme, ngen, UG, OV; heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = allelefreqs_p, pf = ploidy_freq)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[1]]
	p3 = [ploidy_freq(d)[2]]
	p4 = [ploidy_freq(d)[3]]
	
	for n=1:ngen
		d = mating_PnB(d,UG,OV)
		#d = mutate(d) 
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[1])
		push!(p3, ploidy_freq(d)[2])
		push!(p4, ploidy_freq(d)[3])
		
	end
	(het=het, pop=pop, tm=tm, af=af, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen)
end

# ╔═╡ 684280c2-7160-11eb-3353-2f9f783d2ea3
sim_ploidyvar = evolving_deme_ploidyvar(d_p,t,UG,OV)

# ╔═╡ 2f3a10b0-717c-11eb-129b-1d013132c637
ploidy.(sim_ploidyvar.deme.agents)

# ╔═╡ 3394adb0-717b-11eb-1adc-8bdeb0dde24b
begin
	pf2_p = sim_ploidyvar.p2
	pf3_p = sim_ploidyvar.p3
	pf4_p = sim_ploidyvar.p4
	plot(pf2_p, grid=false, color=:blue, label="diploid")
	plot!(pf3_p, grid=false, color=:green, label="triploid")
	plot!(pf4_p, grid=false, color=:red, label="tetraploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ fc62edb0-717b-11eb-2a0c-35af2be1d06c
begin
	popsize_p = map(mean, sim_ploidyvar.pop)
	plot(popsize_p, grid=false, color=:black, label=false, title="Simulated population size")
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ d53330b0-717b-11eb-09dd-c97c372c04f0
begin
	traitmean_ploidy = map(mean, sim_ploidyvar.tm)
	plot(traitmean_ploidy, grid=false, color=:black, label=false, title="Effect of stabilizing selection")
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ fa764302-71f6-11eb-15c6-cfcb5734010b
begin
	Hₒ_ploidyvar = map(mean, sim_ploidyvar.het)
	plot(Hₒ_ploidyvar, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH mixed ploidy with selection")
	plot!(0:sim_ploidyvar.ngen+1, 
		t->expected_heterozygosity(p*(1-p), t, N),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 3147a3b0-4097-11eb-34be-bf35f43bd3e1
md""" ## Mutation To be continued"""

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

# ╔═╡ 926a0e2e-7c1a-11eb-3c33-cdaa406618c4
md""" ## Multiple demes"""

# ╔═╡ 9b9bac30-342f-11eb-120b-9d2e7354f241
"""
    Habitat{D}
A 1-dimensional habitat, i.e. an array of connected demes. This implements
the migration aspects of the population genetic environment.
"""
@with_kw struct Habitat{D,T}
    demes::Vector{D}
    σ ::T = 1/2 #variance of dispersal
    b ::T = 0.1 #steepness of linear gradient
	θ ::T = 12.5 #phenotypic optimum in the center
    Dm::T = 250. #number of demes to initialize
end

# ╔═╡ dc976fe0-7c1d-11eb-2224-33418463cae8
d_p

# ╔═╡ a1bee470-7c1d-11eb-019e-09dcff0c5c73
habi = Habitat(demes=[d_p])

# ╔═╡ 20f9a2ae-7c24-11eb-2a8a-493c2e8b3dd5
(h::Habitat)(demes) = Habitat(h.demes, h.σ, h.θ, h.b, h.dm)

# ╔═╡ ce14f100-4012-11eb-0daf-19617768225e
function random_walk(h::Habitat, p)
	hab = Habitat(demes=[MixedPloidyDeme(agents=(randagent_p(0.5, 0.5, 50, [2], 0, d=2.)), θ=i.θ) for i in h.demes])
    newdemes = similar(h.demes)
    for (i, deme) in enumerate(h.demes)
        for agent in deme.agents
            step = rand() < p ? rand([-1,1]) : 0 
            if step == -1 && i == 1
                step = 0
            elseif step == 1  && i == length(h)
                step = 0
            end
            push!(hab.demes[i+step], agent)
        end
    end
    #Habitat(demes=newdemes, h.σ, h.b, h.θ, h.Dm)
	hab
end

# ╔═╡ d124a200-4012-11eb-3644-2be0a62eab74
function linear_gradient(h::Habitat) 
	#KK = [i*b for i in 0:Dm-1]
	#KK = [θ for i in 0:Dm-1]
	a = -(h.Dm/2)*h.b + h.θ #where Dm is number of demes and b is gradient
	KK = [(i*h.b)+a for i in 0:h.Dm-1]
    return KK
end

# ╔═╡ d98408f0-4012-11eb-3515-bfd13ab42f56
g_lin = linear_gradient(habi)

# ╔═╡ f9474a80-4012-11eb-1b5c-afa826771136
function initiate_habitat(gradient)
	"""The aim should be to initiate a population for nd_s demes on a linear gradient (with slope b) and with "optimal genetic variance" where one half of the genes are adapted, meaning their clines take the form and spacing as assumed for the deterministic model under linkage equilibrium"""
	
	hab = Habitat(demes=[MixedPloidyDeme(agents=(randagent_p(0.5, 0.5, 50, [2], 0, d=2.)), θ=i) for i in gradient])
	for a in randagent_p(0.5, 0.5, 50, [2], 100, d=2.)
	push!(hab.demes[Int(hab.Dm/2)].agents, a)
	end

	return hab
end

# ╔═╡ ff24c180-4012-11eb-3edb-09d780162b72
hab = initiate_habitat(g_lin)

# ╔═╡ eb36ad90-7c27-11eb-0ab3-33969b64dc36
newdemes = similar(hab.demes)

# ╔═╡ 39a60612-7c28-11eb-324a-238466b5c0aa
hab.demes[1]

# ╔═╡ 00f20770-4013-11eb-172e-7d4cfede7a85
function evolving_habitat(h::Habitat{D}, ngen, UG, OV) where D
	for n = 1:ngen
		h = random_walk(h,0.5)
		#h = Gaussian_dispersal(h,σ)
		new_h = Vector{D}(undef, length(h))
		for (i, d) in enumerate(h.demes)
			d = mating_PnB(d,UG,OV)
			#d = mutate(d)
			new_h[i] = d
		end
		h = Habitat(demes=new_h)
	end
	(h=h, ngen=ngen)
end

# ╔═╡ 08ce7910-7c27-11eb-1674-f97c8db9a82d
Base.length(h::Habitat) = length(h.demes)

# ╔═╡ 0832b7f0-4013-11eb-0d7e-d3f77e87aa9b
sim_hab = evolving_habitat(hab, 500, UG, OV)

# ╔═╡ f39e2370-7c2a-11eb-0cd6-3775f3c203ed
pop_sizes = [length(deme) for deme  in sim_hab[1].demes]

# ╔═╡ 6c282e70-7c2c-11eb-2a89-8971cf1b19f0
begin 
ppf1 = [ploidy_freq(deme)[1] for deme  in sim_hab[1].demes]
ppf2 = [ploidy_freq(deme)[2] for deme  in sim_hab[1].demes]
ppf3 = [ploidy_freq(deme)[3] for deme  in sim_hab[1].demes]
end

# ╔═╡ ce005ba0-7c2b-11eb-0037-655fb8d4d398
begin
	K = 100
	σ = 0.5
	b = 0.1
	Vs = 0.5
	rm = 1.06
	Dm = 250
	s = 1
end

# ╔═╡ 3238bff0-7c2b-11eb-1e9e-07866c41d6f4
margin = (sqrt(2)*b*σ)/((2*rm*sqrt(Vs))-b*σ) .>= 0.15.*pop_sizes.*σ*sqrt(s)

# ╔═╡ 0f7d5380-4013-11eb-2a3b-614df9653550
begin
	pop_sizes_p = map(mean, pop_sizes)
	p1 = plot(pop_sizes_p, grid=false, color=:black, label=false, legendfontsize=5)
	hline!([K], label = "K")
	hline!([K*(1-(σ*b)*(1/(2*sqrt(Vs)*rm)))], label = "Expected pop size")
	vline!([Dm/2], label = "Starting deme")
	plot!([margin]*10, label = "Deterministic range margin")
	plot!(ppf1, grid=false, color=:blue, label="Diploid")
	plot!(ppf2, grid=false, color=:green, label="Triploid")
	plot!(ppf3, grid=false, color=:red, label="Tetraploid")
	xlabel!("Space")
	ylabel!("Population size N")
end

# ╔═╡ 0cf49a80-7cd4-11eb-2837-4bdbf5b5b319
begin
	anim_range = @animate for i ∈ 1:10
		sim_habA = evolving_habitat(hab,i,UG,OV)
		#if i != 1
		#sim_habA = evolving_habitat(sim_habA[1],1,1.06,0.5,10^-6,0.50)
		#end
		pop_sizes = [length(deme) for deme  in sim_habA[1].demes]
		pop_sizes_p = map(mean, pop_sizes)
		ppf1 = [ploidy_freq(deme)[1] for deme  in sim_habA[1].demes]
		ppf2 = [ploidy_freq(deme)[2] for deme  in sim_habA[1].demes]
		ppf3 = [ploidy_freq(deme)[3] for deme  in sim_habA[1].demes]
		p1 = plot(pop_sizes_p, grid=false, color=:black, label=false, legendfontsize=5)
		hline!([K], label = "K")
		hline!([K*(1-(σ*b)*(1/(2*sqrt(Vs)*rm)))], label = "Expected pop size")
		vline!([Dm/2], label = "Starting deme")
		#plot!([margin]*10, label = "Deterministic range margin")
		plot!(ppf1, grid=false, color=:blue, label="Diploid")
		plot!(ppf2, grid=false, color=:green, label="Triploid")
		plot!(ppf3, grid=false, color=:red, label="Tetraploid")
		xlabel!("Space")
		ylabel!("Population size N")
	end every 1
	gif(anim_range, "fizzypoly.gif", fps = 3)
end

# ╔═╡ 375696f0-4013-11eb-1779-677909ec6e11
function f_trait_agents(sim_hab)
	trait_agents = []
	cordst = []
	for (i, deme) in enumerate(sim_hab[1].demes)
		for agent in deme.agents
			t = trait(agent)
			p = (i,t)
			push!(cordst,i)
			push!(trait_agents,t)
		end
	end
	trait_agents, cordst
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

# ╔═╡ a7efee70-7fa6-11eb-3811-55a84a7697b7
function f_het_demes(sim_hab)
	het_demes = []
	cordsh = []
	for (i, deme) in enumerate(sim_hab[1].demes)
	if length(deme) != 0
		het = 0.5^2*sum(heterozygosities_p(deme))
		push!(cordsh,i)
		push!(het_demes,het)
	else
		push!(cordsh,i)
		push!(het_demes,0)
		end
	end
	het_demes, cordsh
end

# ╔═╡ a1fc5762-7fa6-11eb-2b9f-799f43e01b8a
begin
	#het_means_p = map(mean, het_demes)
	het_demes, cordsh = f_het_demes(sim_hab)

	p3 = plot(cordsh, het_demes, grid=false, color=:black, label="Vg_mean deme")
	hline!([0.1*0.5*sqrt(0.5)], label = "E(V_G)") #b*σ*sqrt(Vs)
	xlabel!("Space")
	ylabel!("\$V_G\$")
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
# ╟─a444e532-4097-11eb-0a27-89e61edc9802
# ╠═473ff8d0-4186-11eb-25a9-f71f8fdfd1d9
# ╠═3637a430-4199-11eb-020d-951aeb617f11
# ╠═4903f02e-3e76-11eb-3d88-d1566b3a3523
# ╠═4e2a4c30-3e76-11eb-0d8f-1b919ab9d04b
# ╠═59c7ee30-3e76-11eb-07ee-bb3b33a9f6db
# ╠═7f38e42e-3e76-11eb-26b1-9d54c06c0c34
# ╠═ea898c72-40c5-11eb-233f-074fe4f50a3c
# ╠═7d3cfed0-4184-11eb-2b59-29c1afda519d
# ╠═88e341b0-3e76-11eb-05d0-9d9b75e16a6a
# ╠═9eee9130-3e76-11eb-21e2-4516ca929e40
# ╠═a3966e60-3e76-11eb-2e07-3f8537c4a8b4
# ╠═a87f4c30-3e76-11eb-3ad3-33bf4d5d9089
# ╠═90693580-72a0-11eb-2127-e1ca7c4c1cbc
# ╠═5666a9a0-729e-11eb-22ee-cfd529fe980b
# ╠═55d3df82-729e-11eb-23b0-03eeb621d7b2
# ╠═aed0bab0-3e76-11eb-338f-efff9f2811e6
# ╠═b01bb140-3e76-11eb-2716-958dbfc34606
# ╠═d8a5b240-4093-11eb-1b6e-f1129b1ac5b5
# ╟─b5f6d590-717f-11eb-01dc-fb779d0a2410
# ╠═b4ea1900-7099-11eb-3159-453d15bc2f7a
# ╠═b4e7f620-7099-11eb-0954-f3d8ab8d5103
# ╠═a4b40320-7099-11eb-27cb-1f27ef2645ee
# ╠═940db692-70aa-11eb-1276-9d757f1ed1ba
# ╠═fdd49290-7ba1-11eb-3316-45fa9dc62000
# ╠═2260a980-70af-11eb-2688-4d1b7d6f2f8b
# ╠═23b9e2ce-70ad-11eb-1e5d-1f8ad0c849bb
# ╠═14e8d250-7bbd-11eb-0da7-afeb3fc81d83
# ╠═321bbf50-7bad-11eb-0b6d-6f83a6916170
# ╠═737a7620-7bb3-11eb-14c2-83beec620ea1
# ╠═ad2350a0-7bb2-11eb-2517-1f70162a1a54
# ╠═eabf2cb0-7bb0-11eb-3e49-db4d6b9a21e2
# ╠═44dbde9e-7bb1-11eb-0d1b-4b44d5879262
# ╠═73b40ea0-7bb1-11eb-264e-f7120f56fba4
# ╠═1d6cff90-715a-11eb-2648-c9e1455933bd
# ╠═01be8d80-715b-11eb-3b9d-17a328be8383
# ╠═b0ef24ee-715a-11eb-0dd3-a9be75572bdb
# ╠═313d5100-714b-11eb-0f11-db84567a7e36
# ╠═684280c2-7160-11eb-3353-2f9f783d2ea3
# ╠═2f3a10b0-717c-11eb-129b-1d013132c637
# ╠═3394adb0-717b-11eb-1adc-8bdeb0dde24b
# ╠═fc62edb0-717b-11eb-2a0c-35af2be1d06c
# ╠═d53330b0-717b-11eb-09dd-c97c372c04f0
# ╠═fa764302-71f6-11eb-15c6-cfcb5734010b
# ╟─3147a3b0-4097-11eb-34be-bf35f43bd3e1
# ╠═b42552a0-3e76-11eb-3ee2-a5bf1aabcff6
# ╠═619f8c50-408b-11eb-184a-95d63ec67247
# ╠═bdd862b0-3e76-11eb-09fa-6fa07f6d6573
# ╟─926a0e2e-7c1a-11eb-3c33-cdaa406618c4
# ╠═9b9bac30-342f-11eb-120b-9d2e7354f241
# ╠═dc976fe0-7c1d-11eb-2224-33418463cae8
# ╠═a1bee470-7c1d-11eb-019e-09dcff0c5c73
# ╠═20f9a2ae-7c24-11eb-2a8a-493c2e8b3dd5
# ╠═ce14f100-4012-11eb-0daf-19617768225e
# ╠═eb36ad90-7c27-11eb-0ab3-33969b64dc36
# ╠═39a60612-7c28-11eb-324a-238466b5c0aa
# ╠═d124a200-4012-11eb-3644-2be0a62eab74
# ╠═d98408f0-4012-11eb-3515-bfd13ab42f56
# ╠═f9474a80-4012-11eb-1b5c-afa826771136
# ╠═ff24c180-4012-11eb-3edb-09d780162b72
# ╠═00f20770-4013-11eb-172e-7d4cfede7a85
# ╠═08ce7910-7c27-11eb-1674-f97c8db9a82d
# ╠═0832b7f0-4013-11eb-0d7e-d3f77e87aa9b
# ╟─f39e2370-7c2a-11eb-0cd6-3775f3c203ed
# ╟─6c282e70-7c2c-11eb-2a89-8971cf1b19f0
# ╟─ce005ba0-7c2b-11eb-0037-655fb8d4d398
# ╟─3238bff0-7c2b-11eb-1e9e-07866c41d6f4
# ╟─0f7d5380-4013-11eb-2a3b-614df9653550
# ╟─0cf49a80-7cd4-11eb-2837-4bdbf5b5b319
# ╟─375696f0-4013-11eb-1779-677909ec6e11
# ╟─472b52f0-4013-11eb-1e90-4775755333ef
# ╟─3ddde280-4013-11eb-39a3-130d6ec349ff
# ╟─51427c00-4013-11eb-2be7-511a95abdf44
# ╟─a7efee70-7fa6-11eb-3811-55a84a7697b7
# ╟─a1fc5762-7fa6-11eb-2b9f-799f43e01b8a
