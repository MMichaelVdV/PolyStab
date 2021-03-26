# Basic building blocks (structures):

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
mutation). OV: Viability matrix, a symmetric matrix that contains the viability
of offspring for each possible combination of gametes. Ug: Unreduced gamete
formation matrix, a matrix that contains the probability of unreduced gametes
for each level of ploidy in the population.

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
    α ::T     = 0.25
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
(d::MixedPloidyDeme)(agents) = 
    MixedPloidyDeme(agents, d.K, d.θ, d.rm, d.Vs, d.u, d.α, d.μ, d.OV, d.UG)
(d::SimpleDeme)(agents) = SimpleDeme(agents, d.K)
(h::Habitat)(demes) = Habitat(h.demes, h.σ, h.θ, h.b, h.dm)

Base.getindex(a::Agent, i) = a.loci[i]
Base.getindex(a::Agent, i, j) = a.loci[i][j] 
Base.getindex(d::AbstractDeme, i) = d.agents[i]

"""
    randagent(p, α, n, k; d)

Draw a random agent with `n` bi-allelic `k`-ploid loci of allelic effect `α`
with allele frequency `p`
"""
randagent(p, α, n, k; d=1.) = Agent(rand(Bernoulli(p), k, n) * α, d)
randagent(p, α, n, k, N; d=1.) = [randagent(p, α, n, k, d=d) for i=1:N]

# suggested change: k is a vector of weights for ploidy levels 
"""
    randagent_p(p, α, n, k::Vector, N; d)

Sample a bunch of random agents with ploidy levels distributed with weights
`k`. See also `randagent`.
"""
randagent_p(p, α, n, k, N; d=1.) = map(1:N) do x
    randagent(p, α, n, sample(weights(k)), d=d) 
end

# Example:
# A mixed ploidy deme with ~25 diploids and ~25 tetraploids, where α is 0.5 and
# number of loci is 50.  
# d_p = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 0.5, 0., 0.5]))

Base.rand(rng::AbstractRNG, d::AbstractDeme, n) = rand(rng, d.agents, n)
Base.rand(d::AbstractDeme, n::Int) = rand(d.agents, n)

Base.length(a::Agent) = size(a.loci)[2]	#assumes all chromosomes are same length
Base.length(d::AbstractDeme) = length(d.agents)
Base.length(h::Habitat) = length(h.demes)

# I think using `nloci` instead of length leads to more readable
nloci(a::Agent) = size(a.loci)[2]

Base.push!(d::AbstractDeme, a::Agent) = push!(d.agents, a)

Base.sum(a::Agent) = sum(a.loci)

ploidy(a::Agent) = size(a.loci)[1]

emptycopy(d::MixedPloidyDeme{A,T}) where A where T = 
    MixedPloidyDeme(A[], d.K, d.θ, d.rm, d.Vs, d.u, d.μ, d.OV, d.UG)
emptycopy(h::Habitat) = Habitat(emptycopy.(h.demes), h.σ, h.b, h.θ, h.Dm)

expected_heterozygosity(H₀, t, N) = ((1.0-1.0/N)^t)*H₀
	
#Deme level

#Mating:

"""
	viability(a::Agent, b::Agent, d::deme)
"""
function viability(a::Agent, b::Agent, d::AbstractDeme)
	return d.OV[ploidy(a), ploidy(b)]
end

"""
	recombine_poly(a::Agent)

Free recombination between loci in a mixed ploidy population. 
"""
function recombine_poly(a::Agent)
    genome = similar(a.loci)
    # this was incorrect, you actually resample (with replacement) homologs 
    # from the parental genome, and this will not generally be valid. I.e.
    # in this implementation you can have an `Aaaa` individual recombining 
    # to form an `AAAA` genome. Double reduction complicates this, but in
    # the absence of double reduction (i.e. in the case of strict bivalent
    # formation at meiosis) we should do a shuffling, not a resampling.
    for i in 1:nloci(a)
        genome[:,i] = shuffle(a.loci[:,i])
	end
    Agent(loci=genome, d=a.d) 
    # This was 1. * ploidy(a) but I don't get that? 
    # shouldn't it just be the actual `d`?
end

"""
	recombine_poly(d::AbstractDeme{A})

Free recombination between loci in a mixed ploidy population.
"""
function recombine_poly(d::MixedPloidyDeme{A}) where A
	new_agents = A[]
	for agent in d.agents
        push!(new_agents, recombine_poly(agent))
	end
    d(new_agents)
end	
	
# why is this function named unreducaed gamete? It's confusing as it also
# is used for producing reduced gametes...
"""
	unreduced_gamete(a::Agent, d::AbstractDeme)

Unreduced gamete formation in a mixed ploidy population of 2n,3n,4n as it is
implemented at the moment.
"""
function unreduced_gamete(a::Agent, d::AbstractDeme)
	# this samples the ploidy level (1 to 4, potentially) of gametes 
    # why not include an all-zero row for haploid individuals in `d.UG`?
    # (that way you can index by ploidy directly)
	num = sample(1:4, weights(d.UG[ploidy(a)-1,:]))
    idx = sample(1:ploidy(a), num, replace=false)
    return Agent(a.loci[idx, :], a.d)
end	

# XXX: There seems to be some confusion on gamete formation. We should have
# recombination before gamate formation, and random segregation of homologs
# into gametes. I don't think the current code does this properly. Maybe it
# would be instructive to really implement a full meiosis, we might code 
# some shortcuts later for efficiency, but it might be helpful. What I mean
# is implement meiosis algorithmically:
#
# 1. replicate the genome e.g. for a single diploid locus (Aa) -> (AA|aa)
# 2. recombination: (AA|aa) -> (AA|aa) w.p. 1-c
#                           -> (Aa|Aa) w.p. c
# 3. meiosis 1:   (AA|aa) -> (AA) and (aa) 
#                 (Aa|Aa) -> (Aa) and (Aa) 
# 4. meiosis 2:   (Aa) -> (A) (a), (aa) -> (a) (a), (AA) -> (A) (A)
#
# (in the notation above everything in parenthesis is a cell, everything on
# one side of a `|` are sister chromatids)
# Of course for a diploid, this just results in 50% (a) gametes and 50% (A)
# gametes, so nobody cares. For tetraploids this is different, and the frequency
# of each gamete depends on the mechanism of unreduced gamete formation and
# the frequency of double reduction. We may assume, for starters, that there 
# is no double reduction.

"""
	unreduced_gamete(d::AbstractDeme{A})

Unreduced gamete formation in a mixed ploidy population of 2n,3n,4n as it is
implemented at the moment.
"""
function unreduced_gamete(d::MixedPloidyDeme{A}) where A
	new_agents = A[]
	for agent in d.agents
		push!(new_agents, unreduced_gamete(agent,d))
	end
    d(new_agents)
end
	
"""
	mate_p(a::Agent, b::Agent)

Mating in a mixed ploidy deme. Assumes that different cytotypes can be
compatible with a decrease in viability (cfr. OffspringViability matrix) (i.e.
when #individuals with different ploidy hybridize, they generate have a
probability p to generate no viable offspring).  Selfing is allowed without
cost. This influences the dynamics by including hybrid offspring that can
compete for space (and affects the malthusian fitness/density dependence of
selection).
"""	
function mate_p(a::Agent, b::Agent, d::MixedPloidyDeme)
    #gamete formation
    # XXX we should have recombinatoin *before* gamete formation ?!
    ag = unreduced_gamete(recombine_poly(a), d)
    bg = unreduced_gamete(recombine_poly(b), d)
    #combine gametes and assign viability
    via = viability(ag, bg, d) 
    # not sure if returning 0 is the best idea in terms of type stability
    # perhaps you should return a mock-agent (for instance an agent with length
    # zero genome) when there is no viable offspring
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
    d(new_agents)
end 

"""
"""
function number_of_offspring(d::AbstractDeme, a::Agent)
    logw = malthusian_fitness(d, a)
    rand(Poisson(exp(logw)))
end

"""
	mating_PnB(d::AbstractDeme{A})

Mating in a mixed ploidy deme with unreduced gamete formation and partner
choice weighted by fitness (cfr. PnB).
"""
function mating_PnB(d::AbstractDeme{A}) where A
	new_agents =  A[]
	fitnesses = exp.(malthusian_fitness(d))
	for i=1:length(d)
        # so here an individual has all it's offspring with the same partner?
        # perhaps it would be more reasonable to have `noff` parental pairs
        # for B1? (would make more sense for plants at least?). In that case 
        # we'd have something like
        # Bs = sample(d.agents, weights(fitnesses), noff)
        # new_agents = vcat(new_agents, filter(!ismock, map(B2->mate_p(B1, B2), Bs)))
        # where `ismock` checks whether an offspring is a mock agent (for a failed
        # agent, cfr. remark above, or if you keep the `0`, ismock would be x->x==0)
        # see function `suggested` below
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
    d(new_agents)
end 

ismock(x) = x == 0

function suggested(d::AbstractDeme{A}) where A
	new_agents = A[]
	fitnesses = exp.(malthusian_fitness(d))
	for i=1:length(d)
		B1 = d.agents[i]
		noff = number_of_offspring(d, B1)
        Bs = sample(d.agents, weights(fitnesses), noff) 
        offspring = filter(!ismock, map(B2->mate_p(B1, B2), Bs))
        new_agents = vcat(new_agents, offspring)
    end
    d(new_agents)
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

Initialize a vector of phenotypic optima for each deme of the habibat according
to a linear gradient. 
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

Aim should be to initiate a population for nd_s demes on a linear gradient
(with slope b) and with optimal genetic variance where one half of the genes
are adapted, meaning their clines take the form and spacing as assumed for the
deterministic model under linkage equilibrium.
"""
function initiate_habitat(gradient, d::MixedPloidyDeme, p, α, L, N)
    agents = randagent_p(p, α, L, [0., 1., 0., 0.], 0, d = 2.)
    hab = Habitat(demes=[d(agents) for i in gradient])
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

# α = 0.25 #need to incorporate this somewhere
# just put in deme for now?

# XXX this is very inefficient! in most cases there will be no mutations
# so you are loopingfor nothing!
function _mutate(d::AbstractDeme, a::Agent)
	num = ploidy(a)
	loci = zeros(num, length(a))
    newloci = similar(a.loci)
	for i in 1:ploidy(a)
    	for j in 1:length(a)
    		if rand() > d.μ
            	newloci[i,j] = a.loci[i,j]
        	else
            	a.loci[i,j] == 0.0 ? x = d.α : x = 0.0
            	newloci[i,j] = x
        	end
		end
    end
    Agent(newloci, 1. *num)  # here again, what happens with `d` parameter?
end

function mutate(d::AbstractDeme, a::Agent)
    Nsites = nloci(a) * ploidy(a)
    nmutations = rand(Poisson(Nsites * d.μ))
    nmutations == 0 && return a  # or return a copy?
    mutated = sample(1:Nsites, nmutations, replace=false)
    newloci = copy(a.loci)
    for i in mutated
        newloci[i] = newloci[i] == 0. ? d.α : 0.
    end
    return Agent(newloci, a.d)
end
	
function mutate(d::AbstractDeme{A}) where A
    newdeme = Vector{A}(undef, length(d))
    for i=1:length(d)
		newdeme[i] = mutate(d, d.agents[i])
	end
    d(newdeme)
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

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.
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

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.
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

Simulate a habitat with mixed ploidy, malthusian fitness and unreduced gamete
formation.
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
    length(d.agents) == 0 && return [0.]
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
	freq ./ length(d.agents)
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
    length(d.agents) == 0 && return [0.]
    map(p->p*(1-p), allelefreqs_p(d))
end

"""
	ploidy_freq(d::AbstractDeme)

Random function I quickly wrote to check the population size of haploids and
diploids seperately. Should be generalized.
"""
ploidy_freq(d::AbstractDeme) = counts(map(ploidy, d.agents))

var_add(a::Agent,α) = ploidy(a)*α^2*sum(heterozygosities_p(a))

#Genotype -> phenotype maps

"""
    trait(a::Agent)
"""
trait(a::Agent) = sum(a)/a.d  # is this how we define it? shouldn' ploidy come in there?
	
"""
	trait_mean(d::AbstractDeme)
"""
function trait_mean(d::AbstractDeme)
    length(d.agents) == 0 && return 0.
    z = map(trait, d.agents)
    sum(z)/length(d)
end

"""
	f_trait_agents(d::MixedPloidyDeme)
"""
f_trait_agents(d::MixedPloidyDeme) = map(trait, d.agents)
