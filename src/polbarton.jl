#Basic building blocks (structures):

"""
    Agent{T}
Agent of arbitrary ploidy level.
"""
@with_kw struct Agent{T}
    loci::Matrix{T}
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
mutation). OV: Viability matrix, a symmetric matrix that contains the viability of offspring 
for each possible combination of gametes. Ug: Unreduced gamete formation matrix, 
a matrix that contains the probability of unreduced gametes for each level of ploidy in the population.

- `K` : Carrying capacity
- `θ` : Environmental optimum
- `rm`: Mean Malthusian fitness
- `Vs`: Variance of stabilizing selection
- `u` : Unreduced gamete formation rate
- `μ` : Mutation rate
- `OV` : Offspring viability
- `UG` : Unreduced gamete formation
"""
@with_kw struct MixedPloidyDeme{A,T} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 200
    θ ::T     = 12.5
    rm::T     = 1.06
    Vs::T     = 1/2
    u ::T     = 0.01
    μ ::T     = 1e-6
	OV::Matrix{T} = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.]
	UG::Matrix{T} = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.]
end

"""
    Habitat{D}
A 1-dimensional habitat, i.e. an array of connected demes. This implements
the migration aspects of the population genetic environment.
"""
@with_kw struct Habitat{D,T}
    demes::Vector{D}
    σ ::T = sqrt(1/2) #variance of dispersal
    b ::T = 0.1 #steepness of linear gradient
	θ ::T = 12.5 #phenotypic optimum in the center
    Dm::T = 250. #number of demes to initialize
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
randagent_p(p, α, n, k; d=1.) = Agent((rand(Bernoulli(p), k, n) * α), d)
#(rand(Bernoulli(0.5), 10, 5) * 1)
randagent_p(p, α, n, k, N; d=1.) = [randagent_p(p, α, n, rand(k), d=d) for i = 1:N]

#Example:
#A mixed ploidy deme with 25 diploids and 25 tetraploids, where α is 0.5 and number of loci is 50.
#d_p = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 0.5, 50, [2], 25, d=2.),randagent_p(0.5, 0.5, 50, [4], 25, d=4.)))

Base.rand(rng::AbstractRNG, d::AbstractDeme, n) = rand(rng, d.agents, n)
Base.rand(d::AbstractDeme, n::Int) = rand(d.agents, n)

Base.length(a::Agent) = size(a.loci)[2]	#assumes all chromosomes are same length
Base.length(d::AbstractDeme) = length(d.agents)
Base.length(h::Habitat) = length(h.demes)

Base.push!(d::AbstractDeme, a::Agent) = push!(d.agents, a)

Base.sum(a::Agent) = sum(a.loci)

ploidy(a::Agent) = size(a.loci)[1]
ploidy(d::AbstractDeme) = [length(a) for a in d.agents]

emptycopy(d::MixedPloidyDeme{A,T}) where A where T = MixedPloidyDeme(A[], d.K, d.θ, d.rm, d.Vs, d.u, d.μ, d.OV, d.UG)
emptycopy(h::Habitat) = Habitat(emptycopy.(h.demes), h.σ, h.b, h.θ, h.Dm)

expected_heterozygosity(H₀, t, N) = ((1.0-1.0/N)^t)*H₀
	
#Deme level

#Mating:

"""
	viability(a::Agent, b::Agent, d::deme)
"""
function viability(a::Agent, b::Agent, d::AbstractDeme)
	return d.OV[ploidy(a),ploidy(b)]
end

"""
	recombine_poly(a::Agent)
Free recombination between loci in a mixed ploidy population. 
"""
function recombine_poly(a::Agent)
	num = ploidy(a)
	loci = zeros(num, length(a))
	newlocus = similar(a.loci[1,:])
		
	for l in 1:num
		for j in 1:length(a) #this loops over the different loci for each chrosome
			i = rand([x for x in 1:num])
			@inbounds newlocus[j] = a.loci[i,j]
		end
		loci[l,:] += newlocus
	end
	Agent(loci=loci, d= 1. *num)
end

"""
	recombine_poly(d::AbstractDeme{A})
Free recombination between loci in a mixed ploidy population.
"""
function recombine_poly(d::MixedPloidyDeme{A}) where A
	new_agents =  A[]
	for agent in d.agents
		num = ploidy(agent)
		loci = zeros(num, length(agent))
		newlocus = similar(agent.loci[1,:])
		
		for l in 1:num
			for j in 1:length(agent) #this loops over the different loci for each chrosome
				i = rand([x for x in 1:num])
				@inbounds newlocus[j] = agent.loci[i,j]
			end
			loci[l,:] += newlocus
		end
		
		push!(new_agents, Agent(loci=loci,d= 1. *num))
	end
MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ,OV=d.OV,UG=d.UG)
end	
	
"""
	unreduced_gamete(a::Agent, d::AbstractDeme)
Unreduced gamete formation in a mixed ploidy population of 2n,3n,4n as it is implemented at the moment.
"""
function unreduced_gamete(a::Agent, d::AbstractDeme)
	#this samples the ploidy level (1 to 4, potentially) of gametes 
	num = sample([1.,2.,3.,4.], weights(d.UG[ploidy(a)-1,:]))
	loci = zeros(Int(num), length(a))
	b = deepcopy(a.loci)
	#shuffle!(b)
	c = 1
	while c <= num
		i = b[1,:]
		b = b[1:end .!= c,:]
		loci[c,:] += i
		c += 1
	end
	return Agent(loci, 1. *num)
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
	mate_p(a::Agent, b::Agent)
Mating in a mixed ploidy deme. Assumes that different cytotypes can be compatible with a 
decrease in viability (cfr. OffspringViability matrix) (i.e. when #individuals with 
different ploidy hybridize, they generate have a probability p to generate no viable offspring). 
Selfing is allowed without cost. This influences the dynamics by including hybrid offspring that can 
compete for space (and affects the malthusian fitness/density dependence of selection).
"""	
function mate_p(a::Agent, b::Agent, d::MixedPloidyDeme)
    #gamete formation
    ag = recombine_poly(unreduced_gamete(a,d))
    bg = recombine_poly(unreduced_gamete(b,d))
    #combine gametes and assign viability
    via = viability(ag, bg, d) 
    return rand() < via ? Agent([ag.loci ; bg.loci], 1. *(ploidy(ag) + ploidy(bg))) : 0
end

"""
	random_mating_mixedp(d::AbstractDeme{A}) where A
Random mating in a mixed ploidy deme.
"""
function random_mating_mixedp(d::AbstractDeme{A}) where A
    new_agents =  A[]
    for i=1:length(d)
		pair = mate_p(rand(d, 2)...,d)
		if pair != 0
			push!(new_agents, pair)
		end
	end
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ,OV=d.OV,UG=d.UG)
end 

"""
"""
number_of_offspring(d::AbstractDeme,a::Agent) = rand(Poisson(exp(malthusian_fitness(d::AbstractDeme,a::Agent))))

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

#Habitat level

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
function initiate_habitat(gradient,d::MixedPloidyDeme,p::Float64,α::Float64,L::Int64,N::Int64)
	hab = Habitat(demes=[MixedPloidyDeme(agents=(randagent_p(p, α, L, [2], 0, d = 2.)),θ=i,K=d.K,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ,OV=d.OV,UG=d.UG) for i in gradient])
	for a in randagent_p(p, α, L, [2], N, d = 2.)
	push!(hab.demes[Int(hab.Dm/2)].agents, a)
	end

	return hab
end

#Dispersal:

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

#Mutation:

α = 0.25 #need to incorporate this somewhere

function mutate(d::AbstractDeme, a::Agent)
	num = ploidy(a)
	loci = zeros(num, length(a))
    newloci = similar(a.loci)
	for i in 1:ploidy(a)
    	for j in 1:length(a)
    		if rand() > d.μ
            	newloci[i,j] = a.loci[i,j]
        	else
            	a.loci[i,j] == 0.0 ? x = α : x = 0.0
            	newloci[i,j] = x
        	end
		end
    end
    Agent(newloci, 1. *num)
end
	
function mutate(d::AbstractDeme{A}) where A
    newdeme = Vector{A}(undef, length(d))
    for i=1:length(d)
		newdeme[i] = mutate(d, d.agents[i])
	end
MixedPloidyDeme(agents=newdeme,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ,OV=d.OV,UG=d.UG)
end 

#Simulations:

"""
	neutral_evolving_deme(d::AbstractDeme, ngen)
Simulate a single random mating deme with mixed ploidy.
"""
function neutral_evolving_deme(d::MixedPloidyDeme, ngen; heterozygosities_p = heterozygosities_p, allelefreqs_p = allelefreqs_p, trait_mean = trait_mean, pf = ploidy_freq)
	het = [heterozygosities_p(d)]
	af = [allelefreqs_p(d)]
	tm = [trait_mean(d)]
	p2 = [ploidy_freq(d)[1]]
	p3 = [ploidy_freq(d)[2]]
	p4 = [ploidy_freq(d)[3]]
	
	for n=1:ngen
		d = random_mating_mixedp(d)
		#d = unreduced_gamete(d)
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
	evolving_deme_ploidyvar(d::AbstractDeme, ngen)
Simulate a single deme with mixed ploidy, malthusian fitness and unreduced gamete formation.
"""
function evolving_deme_ploidyvar(d::MixedPloidyDeme, ngen; 
	heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = 
	allelefreqs_p, pf = ploidy_freq, fta = f_trait_agents)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[1]]
	p3 = [ploidy_freq(d)[2]]
	p4 = [ploidy_freq(d)[3]]
	fta = [f_trait_agents(d)]
	
	for n=1:ngen
		d = mating_PnB(d)
		d = mutate(d) 
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[1])
		push!(p3, ploidy_freq(d)[2])
		push!(p4, ploidy_freq(d)[3])
		push!(fta, f_trait_agents(d))
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, het=het,tm=tm, af=af, fta=fta) 
end
	
"""
	evolving_deme_UG(d::AbstractDeme, ngen)
Simulate a single deme with mixed ploidy, malthusian fitness and unreduced gamete formation.
"""
function evolving_deme_UG(d::MixedPloidyDeme, ngen; pf = ploidy_freq)

	pop = [length(d)]
	p2 = [ploidy_freq(d)[1]]
	p3 = [ploidy_freq(d)[2]]
	p4 = [ploidy_freq(d)[3]]
	
	for n=1:ngen
		d = mating_PnB(d)
		d = mutate(d) 
		push!(pop, length(d))
		push!(p2, ploidy_freq(d)[1])
		push!(p3, ploidy_freq(d)[2])
		push!(p4, ploidy_freq(d)[3])
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen) 
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
			d = mutate(d)
			new_h[i] = d
		end
		h = Habitat(demes=new_h)
		push!(data, h)
	end
	(h=h, ngen=ngen, data=data)
end

#Utility and plotting

"""
	allelefreqs_p(d::AbstractDeme)
"""
function allelefreqs_p(d::AbstractDeme)
	if length(d.agents) > 0
	freq = Vector{Float64}(undef,length(d.agents[1]))
	for loc in 1:length(freq)
		s = 0
    	for ag in d.agents
			for chr in 1:ploidy(ag)
        		if ag.loci[chr, loc] != 0
            		s += 1/ploidy(ag)
				end
			end
        end
		f = s
		freq[loc] = f
	end
	freq./length(d.agents)
	else [0.]
	end
end	

"""
	allelefreqs_p(a::Agent)
"""
function allelefreqs_p(a::Agent)
	freq = Vector{Float64}(undef,length(a))
	for loc in 1:length(freq)
		s = 0.
    	for chr in 1:ploidy(a)
        	if a.loci[chr, loc] != 0.
            	s += 1. /ploidy(a)
			end
		end
		freq[loc] = s
	end
	freq
end	
	
"""
	heterozygosities_p(d::AbstractDeme)
"""
function heterozygosities_p(d::AbstractDeme)
	if length(d.agents) > 0
	freqs=allelefreqs_p(d)
	map(p->p*(1-p), freqs)
	else [0.]
	end
end

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

var_add(a::Agent,α) = ploidy(a)*α^2*sum(heterozygosities_p(a))

#Genotype -> phenotype maps

"""
"""
trait(a::Agent) = sum(a)/a.d
	
"""
	trait_mean(d::AbstractDeme)
"""
function trait_mean(d::AbstractDeme)
	if length(d.agents) > 0
	z = Float64[]
	for agent in d.agents
		push!(z,trait(agent))
	end
	sum(z)/length(d)
	else 0.
	end
end

"""
	f_trait_agents(d::MixedPloidyDeme)
"""
function f_trait_agents(d::MixedPloidyDeme)
	trait_agents = Float64[]
	for agent in d.agents
		t = trait(agent)
		push!(trait_agents,t)
	end
	trait_agents
end
