### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 58c46a2e-85ee-11eb-28ed-f927eac6cd9b
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes

# ╔═╡ 4d51f440-8659-11eb-010c-ab5401fcc485
md""" ###### This notebook contains the code for generating GIFs for mixed ploidy migration along a 1D habitat."""

# ╔═╡ dac9ff40-8657-11eb-380c-23337b4e14a2
md""" ### Functions
"""

# ╔═╡ a6e75fb0-8657-11eb-19c0-73a0c663c362
begin

#Basic building blocks (structures):

"""
    Agent{T}
Agent of arbitrary ploidy level.
"""
@with_kw struct Agent{T,N}
    loci::Array{Array{T,1},N}
    d::T = 1. # allelic effect scaler for different ploidy levels
end

"""
    AbstractDeme{A}
"""
abstract type AbstractDeme{A} 
end

"""
    SimpleDeme{A}
A single random-mating, single ploidy level deme. 
- `K` : Carrying capacity
"""
@with_kw struct SimpleDeme{A} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 15 
end
	
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
- `OV` : Offspring viability
- `UG` : Unreduced gamete formation
"""
@with_kw struct MixedPloidyDeme{A,T,N} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 100
    θ ::T     = 12.5
    rm::T     = 1.06
    Vs::T     = 1/2
    u ::T     = 0.01
    μ ::T     = 1e-6
	OV::Array{Array{T,1},N} = [[1.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]]
	UG::Array{Array{T,1},N} = [[1.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]]
end

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

"""
    OffspringViability{T,N}
Viability matrix, a symmetric matrix that contains the viability of offspring for each possible combination of gametes. 
(should probably become integrated within other structure,i.e. at either agent or deme level)
"""	
struct OffspringViability{T,N}
    viability::Array{Array{T,1},N}
end

"""
    UnreducedGamete{T,N} 
Unreduced gamete formation matrix, a matrix that contains the probability of unreduced gametes for each level of ploidy in the population.
(should probably become integrated within other structure,i.e. at either agent or deme level)
"""	
struct UnreducedGamete{T,N} 
    prob::Array{Array{T,1},N}
end

#Some useful short functions:

(d::MixedPloidyDeme)(agents) = MixedPloidyDeme(agents, d.K, d.θ, d.rm, d.Vs, d.u, d.μ, d.OV, d.UG)
(d::SimpleDeme)(agents) = SimpleDeme(agents, d.K)
(h::Habitat)(demes) = Habitat(h.demes, h.σ, h.θ, h.b, h.dm)

Base.getindex(a::Agent, i) = a.loci[i]
Base.getindex(a::Agent, i, j) = a.loci[i][j] 
Base.getindex(d::AbstractDeme, i) = d.agents[i]

randagent(p, α, n; d=1.) = Agent([rand(Bernoulli(p), n) * α], d)
randagent(p, α, n, N; d=1.) = [randagent(p, α, n, d=d) for i=1:N]

#'k' is a vector of different ploidy levels to intiatiate a mixed ploidy population
randagent_p(p, α, n, k; d=1.) = Agent([(rand(Bernoulli(p), n) * α) for i=1:k], d)
randagent_p(p, α, n, k, N; d=1.) = [randagent_p(p, α, n, rand(k), d=d) for i=1:N]

#Example:
#A mixed ploidy deme with 25 diploids and 25 tetraploids, where α is 0.5 and number of loci is 50.
#d_p = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 0.5, 50, [2], 25, d=2.),randagent_p(0.5, 0.5, 50, [4], 25, d=4.)))

Base.rand(rng::AbstractRNG, d::AbstractDeme, n) = rand(rng, d.agents, n)
Base.rand(d::AbstractDeme, n::Int) = rand(d.agents, n)

Base.length(a::Agent) = length(a.loci[1]) #should be ok if all chromosomes are same length
Base.length(d::AbstractDeme) = length(d.agents)
Base.length(h::Habitat) = length(h.demes)

Base.push!(d::AbstractDeme, a::Agent) = push!(d.agents, a)

Base.sum(a::Agent) = sum(sum(a.loci))

ploidy(a::Agent) = length(a.loci)
ploidy(d::AbstractDeme) = [length(a.loci) for a in d.agents]

emptycopy(d::MixedPloidyDeme{A,T}) where A where T = MixedPloidyDeme(A[], d.K, d.θ, d.rm, d.Vs, d.u, d.μ, d.OV, d.UG)
emptycopy(h::Habitat) = Habitat(emptycopy.(h.demes), h.σ, h.b, h.θ, h.Dm)

#Functions on deme level

#Functions for mating:


"""
	mate_p(a::Agent, b::Agent)
Mating in a mixed ploidy deme. Assumes that different cytotypes can be compatible with a 
decrease in viability (cfr. OffspringViability matrix) (i.e. when #individuals with 
different ploidy hybridize, they generate have a probability p to generate no viable offspring). 
Selfing is allowed without cost. This influences the dynamics by including hybrid offspring that can 
compete for space (and affects the malthusian fitness/density dependence of selection).
"""	
function mate_p(a::Agent, b::Agent, d::AbstractDeme)
	#gamete formation
	ag = recombine_poly(unreduced_gamete(a,d))
	bg = recombine_poly(unreduced_gamete(b,d))
	gam_a = ag.loci	
	gam_b = bg.loci
	#combine gametes and assign viability
	via = viability(ag,bg,d)
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

"""
	mating_PnB(d::AbstractDeme{A})
Mating in a mixed ploidy deme with unreduced gamete formation and partner choice weighted by fitness (cfr. PnB).
"""
function mating_PnB(d::AbstractDeme{A}) where A
	new_agents =  A[]
	fitnesses = exp.(malthusian_fitness(d))
	for i=1:length(d)
		B1 = d.agents[i]
		noff = number_of_offspring(d,B1)
		B2 = sample(d.agents, weights(fitnesses))
		#B2 = rand(d.agents)
		#child = mate(B1,B2)
		m = mate_p(B1,B2,d)
		if m != 0
			for c in 1:noff 
				push!(new_agents, m)
			end
		end
	end
MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ,OV=d.OV,UG=d.UG)
end 

"""
	random_mating_mixedp(d::AbstractDeme{A}) where A
Random mating in a mixed ploidy deme.
"""
function random_mating_mixedp(d::AbstractDeme{A}) where A
    new_agents =  A[]
    for i=1:length(d)
		pair = mate_p(rand(d, 2)...)
		if pair != 0
			push!(new_agents, pair)
		end
	end
MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ,OV=d.OV,UG=d.UG)
end 
	
"""
"""
trait(a::Agent) = sum(a)/a.d
	
"""
	trait_mean(d::AbstractDeme)
"""
function trait_mean(d::AbstractDeme)
	z = Float64[]
	for agent in d.agents
		push!(z,trait(agent))
	end
	sum(z)/length(d)
end

"""
	unreduced_gamete(a::Agent, d::AbstractDeme)
Unreduced gamete formation in a mixed ploidy population of 2n,3n,4n as it is implemented at the moment.
"""
function unreduced_gamete(a::Agent, d::AbstractDeme)
	#this samples the ploidy level (1 to 4, potentially) of gametes 
	num = sample([1.,2.,3.,4.],weights(d.UG[ploidy(a)-1]))
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

"""
	unreduced_gamete(d::AbstractDeme{A})
Unreduced gamete formation in a mixed ploidy population of 2n,3n,4n as it is implemented at the moment.
"""
function unreduced_gamete(d::MixedPloidyDeme{A}) where A
	new_agents =  A[]
	for agent in d.agents
		push!(new_agents, unreduced_gamete(agent,d))
	end
MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ,OV=d.OV,UG=d.UG)
end	


"""
	malthusian_fitness(d::AbstractDeme,a::Agent)
"""
function malthusian_fitness(d::AbstractDeme,a::Agent)
    N = length(d)
    z = trait(a)
    return d.rm*(1-(N/d.K))-((z-d.θ)^2)/(2*d.Vs)
end 

"""
	malthusian_fitness(d::AbstractDeme)
"""
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


number_of_offspring(d::AbstractDeme,a::Agent) = rand(Poisson(exp(malthusian_fitness(d::AbstractDeme,a::Agent))))

"""
	viability(a::Agent, b::Agent, d::deme)
"""
function viability(a::Agent, b::Agent, d::AbstractDeme)
	return d.OV[ploidy(a)][ploidy(b)]
end

"""
	recombine_poly(a::Agent)
Free recombination between loci in a mixed ploidy population. 
"""
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

"""
	recombine_poly(d::AbstractDeme{A})
Free recombination between loci in a mixed ploidy population.
"""
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
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ,OV=d.OV,UG=d.UG)
end	

#Functions on habitat level

"""
	linear_gradient(h::Habitat)
Initialize a vector of phenotypic optima for each deme of the habibat according to a linear gradient. 
"""
function linear_gradient(h::Habitat) 
	#KK = [i*b for i in 0:Dm-1]
	#KK = [θ for i in 0:Dm-1]
	a = -(h.Dm/2)*h.b + h.θ #where Dm is number of demes and b is gradient
	KK = [(i*h.b)+a for i in 0:h.Dm-1]
    return KK
end

"""
	initiate_habitat(gradient::Vector)
Aim should be to initiate a population for nd_s demes on a linear gradient (with slope b) and with 
optimal genetic variance where one half of the genes are adapted, meaning their clines take the 
form and spacing as assumed for the deterministic model under linkage equilibrium.
"""
function initiate_habitat(gradient,d::MixedPloidyDeme)
	hab = Habitat(demes=[MixedPloidyDeme(agents=(randagent_p(0.5, 0.5, 50, [2], 0, d=2.)),θ=i,K=d.K,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ,OV=d.OV,UG=d.UG) for i in gradient])
	for a in randagent_p(0.5, 0.5, 50, [2], 100, d=2.)
	push!(hab.demes[Int(hab.Dm/2)].agents, a)
	end

	return hab
end

#Functions for dispersal:

"""
	random_walk(h::Habitat, p)
"""
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

"""
	Gaussian_dispersal(h::Habitat,σ)
"""
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
            push!(new_h.demes[i+step].agents, agent)
        end
    end
    new_h
end

"""
	Cauchy_dispersal(h::Habitat,σ)
"""
function Cauchy_dispersal(h::Habitat,σ)
    new_h = emptycopy(h)
	dist = Cauchy(0,σ)
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

#Functions used for simulations:

"""
	neutral_evolving_deme(d::AbstractDeme, ngen)
Simulate a single random mating deme with mixed ploidy.
"""
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

"""
	evolving_deme_popvar(d::AbstractDeme, ngen)
Simulate a single deme with mixed ploidy and malthusian fitness.
"""
function evolving_deme_popvar(d::AbstractDeme, ngen; heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = allelefreqs_p, pf = ploidy_freq)
	#het = [heterozygosities_p(d)]
	pop = [length(d)]
	#tm = [trait_mean(d)]
	#af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[1]]
	p3 = [ploidy_freq(d)[2]]
	p4 = [ploidy_freq(d)[3]]
	
	for n=1:ngen
		d = mating_PnB(d)
		d = unreduced_gamete(d)
		#d = mutate(d) #is this defined on level of deme or agent ?
		#push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		#push!(tm, trait_mean(d))
		#push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[1])
		push!(p3, ploidy_freq(d)[2])
		push!(p4, ploidy_freq(d)[3])
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen)#het, tm=tm, af=af, 
end

"""
	evolving_deme_ploidyvar(d::AbstractDeme, ngen, UG, OV)
Simulate a single deme with mixed ploidy, malthusian fitness and unreduced gamete formation.
"""
function evolving_deme_ploidyvar(d::AbstractDeme, ngen, UG, OV; heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = allelefreqs_p, pf = ploidy_freq)
	#het = [heterozygosities_p(d)]
	pop = [length(d)]
	#tm = [trait_mean(d)]
	#af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[1]]
	p3 = [ploidy_freq(d)[2]]
	p4 = [ploidy_freq(d)[3]]
	
	for n=1:ngen
		d = mating_PnB(d,UG,OV)
		#d = mutate(d) 
		#push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		#push!(tm, trait_mean(d))
		#push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[1])
		push!(p3, ploidy_freq(d)[2])
		push!(p4, ploidy_freq(d)[3])
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen) #het=het,tm=tm, af=af, 
end

"""
	evolving_deme_ploidyvar(d::AbstractDeme, ngen)
Simulate a single deme with mixed ploidy, malthusian fitness and unreduced gamete formation.
"""
function evolving_deme_ploidyvar(d::MixedPloidyDeme, ngen; heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = allelefreqs_p, pf = ploidy_freq)
	#het = [heterozygosities_p(d)]
	pop = [length(d)]
	#tm = [trait_mean(d)]
	#af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[1]]
	p3 = [ploidy_freq(d)[2]]
	p4 = [ploidy_freq(d)[3]]
	
	for n=1:ngen
		d = mating_PnB(d)
		#d = mutate(d) 
		#push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		#push!(tm, trait_mean(d))
		#push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[1])
		push!(p3, ploidy_freq(d)[2])
		push!(p4, ploidy_freq(d)[3])
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen) #het=het,tm=tm, af=af, 
end

"""
	evolving_habitat(h::Habitat{D}, ngen, UG, OV)
Simulate a habitat with mixed ploidy, malthusian fitness and unreduced gamete formation.
"""
function evolving_habitat(h::Habitat{D}, ngen, UG, OV) where D
	for n = 1:ngen
		h = random_walk(h,0.5)
		#ih = Gaussian_dispersal(h,σ)
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

"""
	evolving_habitat(h::Habitat{D}, ngen)
Simulate a habitat with mixed ploidy, malthusian fitness and unreduced gamete formation.
"""
function evolving_habitat(h::Habitat{D}, ngen) where D
	data = []
	for n = 1:ngen
		h = random_walk(h,0.5)
		#ih = Gaussian_dispersal(h,σ)
		new_h = Vector{D}(undef, length(h))
		for (i, d) in enumerate(h.demes)
			d = mating_PnB(d)
			#d = mutate(d)
			new_h[i] = d
		end
		h = Habitat(demes=new_h)
		push!(data, h)
	end
	(h=h, ngen=ngen, data=data)
end

#Utility

"""
	allelefreqs(d::AbstractDeme)
"""
function allelefreqs(d::AbstractDeme)
    # assuming 0 is absence of allele
    f = j->mapreduce(i->d[i][j] != 0., +, 1:length(d))/length(d)
	frequencies = map(f, 1:length(d[1]))
end

"""
	allelefreqs_p(d::AbstractDeme)
"""
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

"""
	heterozygosities(d::AbstractDeme) 
"""
function heterozygosities(d::AbstractDeme) 
	freqs=allelefreq(d)
    map(p->p*(1-p), freqs)
end

"""
	heterozygosities_p(d::AbstractDeme)
"""
function heterozygosities_p(d::AbstractDeme)
	freqs=allelefreqs_p(d)
    map(p->p*(1-p), freqs)
end

expected_heterozygosity(H₀, t, N) = ((1.0-1.0/N)^t)*H₀

"""
	ploidy_freq(d::AbstractDeme)
Random function I quickly wrote to check the population size of haploids and diploids seperately. Should be generalized.
"""
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


end

# ╔═╡ 2664cbb0-85ef-11eb-253a-d9f02c63f5f9
d_p = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 0.5, 50, [2], 100, d=2.),randagent_p(0.5, 0.5, 50, [4], 0, d=4.)), OV = [[1.,0.1,0.,0.],[0.1,1.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]], UG = [[0.9,0.1,0.,0.],[0.,0.,0.,0.],[0.,0.9,0.,0.1]] )

# ╔═╡ 9e9e62c0-865e-11eb-134c-ef481e41ab27
d_p

# ╔═╡ 9e2fdbc0-8659-11eb-1d53-450c776a5ad0
sim_ploidyvar1 = evolving_deme_ploidyvar(d_p,20)

# ╔═╡ c4e53bce-865d-11eb-304b-25da21c3de6f
begin
	pf2_p1 = sim_ploidyvar1.p2
	pf3_p1 = sim_ploidyvar1.p3
	pf4_p1 = sim_ploidyvar1.p4
	plot(pf2_p1, grid=false, color=:blue, label="diploid")
	plot!(pf3_p1, grid=false, color=:green, label="triploid")
	plot!(pf4_p1, grid=false, color=:red, label="tetraploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end


# ╔═╡ e162ef50-865d-11eb-1a6d-27208c981720
habi = Habitat(demes=[d_p])

# ╔═╡ c3df6060-8664-11eb-2912-917d88e0a358
g_lin = linear_gradient(habi)

# ╔═╡ d0625950-8664-11eb-0ef8-696f5caffc41
hab = initiate_habitat(g_lin, d_p)

# ╔═╡ afa35910-8819-11eb-1aa9-9f350178f789
hab.demes[1].UG

# ╔═╡ 0c3b64d0-8665-11eb-2599-c175731347d3
begin
	K = d_p.K
	σ = hab.σ
	b = hab.b
	Vs = d_p.Vs
	rm = d_p.rm
	Dm = hab.Dm
	s = 1
end

# ╔═╡ e2f8cd60-8664-11eb-02c0-772862e46f2b
sim_hab = evolving_habitat(hab, 100)

# ╔═╡ 962bb660-8668-11eb-0b06-13560951b56b
sim_hab.data

# ╔═╡ f31ba960-8664-11eb-1ea6-f70cc3b57243
pop_sizes = [length(deme) for deme  in sim_hab[1].demes]

# ╔═╡ 111ffce0-8665-11eb-2dcd-75143b7d3432
margin = (sqrt(2)*b*σ)/((2*rm*sqrt(Vs))-b*σ) .>= 0.15.*pop_sizes.*σ*sqrt(s)

# ╔═╡ e7d6d8f0-8668-11eb-32a7-1ffb7fb68ab0
sim_hab.data[1]

# ╔═╡ 07e21cd0-8665-11eb-20a7-9f425dbe0c96
begin 
ppf1 = [ploidy_freq(deme)[1] for deme  in sim_hab[1].demes]
ppf2 = [ploidy_freq(deme)[2] for deme  in sim_hab[1].demes]
ppf3 = [ploidy_freq(deme)[3] for deme  in sim_hab[1].demes]
end

# ╔═╡ 196aca60-8665-11eb-3f6e-63cb1fc323b5
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

# ╔═╡ 2ab891d0-8665-11eb-197e-9b8e15dfe3a2
function f_trait_agents(h::Habitat)
	trait_agents = []
	cordst = []
	for (i, deme) in enumerate(h.demes)
		for agent in deme.agents
			t = trait(agent)
			p = (i,t)
			push!(cordst,i)
			push!(trait_agents,t)
		end
	end
	trait_agents, cordst
end

# ╔═╡ 2c457360-8665-11eb-1dc4-cfe4eec3e5d5
trait_means = [trait_mean(deme) for deme in sim_hab[1].demes]

# ╔═╡ 5bf91c00-8675-11eb-1ad3-9ffa80479b47
function f_het_demes(h::Habitat)
	het_demes = []
	cordsh = []
	for (i, deme) in enumerate(h.demes)
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

# ╔═╡ 47434ea0-881b-11eb-337a-411a3d601eb3
begin
	anim_range = @animate for i ∈ 1:100
		sim_habA = sim_hab.data[i]
		#if i != 1
		#sim_habA = evolving_habitat(sim_habA[1],1,1.06,0.5,10^-6,0.50)
		#end
		pop_sizes = [length(deme) for deme  in sim_habA.demes]
		pop_sizes_p = map(mean, pop_sizes)
		
		ppf1 = [ploidy_freq(deme)[1] for deme  in sim_habA.demes]
		ppf2 = [ploidy_freq(deme)[2] for deme  in sim_habA.demes]
		ppf3 = [ploidy_freq(deme)[3] for deme  in sim_habA.demes]
		het_demes, cordsh = f_het_demes(sim_habA)
		
		plot(pop_sizes_p, grid=false, color=:black, label="Tetraploid")
		#hline!([K], label = "K")
		hline!([K*(1-(σ*b)*(1/(2*sqrt(Vs)*rm)))], label = "Expected pop size")
		#vline!([Dm/2], label = "Starting deme")
		plot!([margin]*5, color=:yellow, label = "Deterministic range margin")
		plot!(ppf1, grid=false, color=:blue, label="Diploid")
		plot!(ppf2, grid=false, color=:green, label="Triploid")
		plot!(ppf3, grid=false, color=:red, label="Tetraploid")
		
		
		xlabel!("Space")
		ylabel!("Population size N")
	end every 1
	gif(anim_range, "fizzypop.gif", fps = 3)
end

# ╔═╡ 737a56be-8682-11eb-3ecb-7d1fc90b66c2
begin
	anim_range_Vg = @animate for i ∈ 1:100
		sim_habA = sim_hab.data[i]
		#if i != 1
		#sim_habA = evolving_habitat(sim_habA[1],1,1.06,0.5,10^-6,0.50)
		#end
		pop_sizes = [length(deme) for deme  in sim_habA.demes]
		pop_sizes_p = map(mean, pop_sizes)
		
		#./pop_sizes
		
		ppf1 = [ploidy_freq(deme)[1] for deme  in sim_habA.demes]/100
		ppf2 = [ploidy_freq(deme)[2] for deme  in sim_habA.demes]/100
		ppf3 = [ploidy_freq(deme)[3] for deme  in sim_habA.demes]/100
		het_demes, cordsh = f_het_demes(sim_habA)
		
		p1 = plot(cordsh, het_demes, grid=false, color=:black, label="Vg_mean deme")
		vline!([Dm/2], label = "Starting deme")
		#plot!([margin]*1, label = "Deterministic range margin")
		plot!(ppf1, grid=false, color=:blue, label="Diploid")
		plot!(ppf2, grid=false, color=:green, label="Triploid")
		plot!(ppf3, grid=false, color=:red, label="Tetraploid")
		
		
		xlabel!("Space")
		ylabel!("Population size N")
	end every 1
	gif(anim_range_Vg, "fizzyVg.gif", fps = 3)
end

# ╔═╡ 409372f0-869b-11eb-3f87-ff7038edde33
begin
	anim_range_trait = @animate for i ∈ 1:100
		sim_habA = sim_hab.data[i]
		trait_means = [trait_mean(deme) for deme in sim_habA.demes]
	    trait_means_p = map(mean, trait_means)
		trait_agents, cordst = f_trait_agents(sim_habA)
		
		plot(g_lin, grid=false, color=:blue, label="Z optimum", linestyle=:dash)
		plot!(cordst,trait_agents, label="Z agents")
		plot!(trait_means_p, grid=false, color=:black, label="Z_mean deme")
		xlabel!("Space")
		ylabel!("Trait Z")
	end every 1
	gif(anim_range_trait, "fizzy.gif", fps = 3)
end

# ╔═╡ Cell order:
# ╟─4d51f440-8659-11eb-010c-ab5401fcc485
# ╠═2664cbb0-85ef-11eb-253a-d9f02c63f5f9
# ╠═9e9e62c0-865e-11eb-134c-ef481e41ab27
# ╠═9e2fdbc0-8659-11eb-1d53-450c776a5ad0
# ╠═c4e53bce-865d-11eb-304b-25da21c3de6f
# ╠═e162ef50-865d-11eb-1a6d-27208c981720
# ╠═c3df6060-8664-11eb-2912-917d88e0a358
# ╠═d0625950-8664-11eb-0ef8-696f5caffc41
# ╠═afa35910-8819-11eb-1aa9-9f350178f789
# ╠═e2f8cd60-8664-11eb-02c0-772862e46f2b
# ╠═962bb660-8668-11eb-0b06-13560951b56b
# ╠═0c3b64d0-8665-11eb-2599-c175731347d3
# ╠═f31ba960-8664-11eb-1ea6-f70cc3b57243
# ╠═07e21cd0-8665-11eb-20a7-9f425dbe0c96
# ╠═111ffce0-8665-11eb-2dcd-75143b7d3432
# ╠═2ab891d0-8665-11eb-197e-9b8e15dfe3a2
# ╠═2c457360-8665-11eb-1dc4-cfe4eec3e5d5
# ╠═196aca60-8665-11eb-3f6e-63cb1fc323b5
# ╠═e7d6d8f0-8668-11eb-32a7-1ffb7fb68ab0
# ╠═5bf91c00-8675-11eb-1ad3-9ffa80479b47
# ╠═47434ea0-881b-11eb-337a-411a3d601eb3
# ╠═737a56be-8682-11eb-3ecb-7d1fc90b66c2
# ╠═409372f0-869b-11eb-3f87-ff7038edde33
# ╟─dac9ff40-8657-11eb-380c-23337b4e14a2
# ╠═58c46a2e-85ee-11eb-28ed-f927eac6cd9b
# ╠═a6e75fb0-8657-11eb-19c0-73a0c663c362
