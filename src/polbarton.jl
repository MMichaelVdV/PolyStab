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
	GametePool{T}
Gamete container.
"""
@with_kw struct GametePool{T}
	gametes::Matrix{T}
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
- `rm`: Growth rate per capita
- `Vs`: Variance of stabilizing selection
- `α` : Allelic effect size
- `μ` : Mutation rate
- `OV` : Offspring viability
- `UG` : Unreduced gamete probabilities
"""
@with_kw struct MixedPloidyDeme{A,T} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 200
    θ ::T     = 12.5
    rm::T     = 1.06
    Vs::T     = 1/2
    α ::T     = 0.25
    μ ::T     = 1e-6
	OV::Matrix{T} = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.]
	UG::Matrix{T} = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.]
end

"""
    IslandDeme{A,T}

A single random-mating, mixed-ploidy level deme, most of the 'population
genetic' environment should be implemented at this level (drift, selection,
mutation). OV: Viability matrix, a symmetric matrix that contains the viability
of offspring for each possible combination of gametes. Ug: Unreduced gamete
formation matrix, a matrix that contains the probability of unreduced gametes
for each level of ploidy in the population.

- `K` : Carrying capacity
- `θ` : Environmental optimum
- `β` : Selection gradient
- `m` : migration rate
- `α` : allelic effect size
- `μ` : Mutation rate
- `OV` : Offspring viability
- `UG` : Unreduced gamete probabilities
"""
@with_kw struct IslandDeme{A,T} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 200
    θ ::T     = 12.5
    β ::T     = 0.1
	m ::T     = 0.1
    α ::T     = 0.25
    μ ::T     = 1e-6
	OV::Matrix{T} = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.]
	UG::Matrix{T} = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.]
end

"""
    Habitat{D}

A 1-dimensional habitat as an array of connected demes.
- `σ` : Variance of dispersal
- `b` : Steepness of linear gradient
- `θ` : Phenotypic optimum in the central deme of the habitat
- `Dm` : Number of demes
"""
@with_kw struct Habitat{D,T}
    demes::Vector{D}
    σ ::T = sqrt(1/2) 
    b ::T = 0.1 
	θ ::T = 12.5 
    Dm::T = 250. 
end

#Some useful short functions:
(a::Agent)(loci)=Agent(loci, a.d)
(d::MixedPloidyDeme)(agents) = 
    MixedPloidyDeme(agents, d.K, d.θ, d.rm, d.Vs, d.α, d.μ, d.OV, d.UG)
(d::IslandDeme)(agents) = 
    IslandDeme(agents, d.K, d.θ, d.β, d.m, d.α, d.μ, d.OV, d.UG)
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
# d_mp = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 0.5, 0., 0.5],50))

Base.rand(rng::AbstractRNG, d::AbstractDeme, n) = rand(rng, d.agents, n)
Base.rand(d::AbstractDeme, n::Int) = rand(d.agents, n)

Base.length(a::Agent) = size(a.loci)[2]	#assumes all chromosomes are same length
Base.length(d::AbstractDeme) = length(d.agents)
Base.length(h::Habitat) = length(h.demes)

# I think using `nloci` instead of length leads to more readable code
nloci(a::Agent) = size(a.loci)[2]
popsize(d::AbstractDeme) = length(d.agents)

Base.push!(d::AbstractDeme, a::Agent) = push!(d.agents, a)

Base.sum(a::Agent) = sum(a.loci)

ploidy(a::Agent) = size(a.loci)[1]

emptycopy(d::MixedPloidyDeme{A,T}) where A where T = 
    MixedPloidyDeme(A[], d.K, d.θ, d.rm, d.Vs, d.α, d.μ, d.OV, d.UG)
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
	recombination(a::Agent)

Free recombination between loci in a mixed ploidy population. 
"""
function recombination(a::Agent)
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
end

"""
	recombination(a::Agent,c)

Recombination between loci with recombination rate `c`  in a mixed ploidy population. 
"""
function recombination(a::Agent,c)
    genome = similar(a.loci)
    for i in 1:nloci(a)
		genome[:,i] = rand() < c ? shuffle(a.loci[:,i]) : a.loci[:,i]
    end
    Agent(loci=genome, d=a.d) 
end

"""
	recombination(d::AbstractDeme{A})

Free recombination between loci in a mixed ploidy population.
"""
function recombination(d::AbstractDeme{A}) where A
	new_agents = A[]
	for agent in d.agents
        push!(new_agents, recombination(agent))
	end
    d(new_agents)
end	
	
"""
	gametogenesis(a::Agent, d::AbstractDeme)

Gamete formation in a mixed ploidy population (1n to 4n).
"""
function gametogenesis(a::Agent, d::AbstractDeme)
	num = sample(1:4, weights(d.UG[ploidy(a),:]))
    idx = sample(1:ploidy(a), num, replace=false)
    return Agent(a.loci[idx, :], a.d)
end	

"""
	meiosis(a::Agent,c)

Meiosis in a diploid.
- `c` : Recombination rate
"""
function meiosis(a::Agent,c)
	# 1. replicate the genome e.g. for a single diploid locus (Aa) -> (AA|aa)
	gametes = [a.loci  a.loci]
	a = a(gametes)
	#2. recombination: (AA|aa) -> (AA|aa) w.p. 1-c
	#                           -> (Aa|Aa) w.p. c
	a = recombination(a,c)
	# 3. meiosis 1:   (AA|aa) -> (AA) and (aa) 
	#                 (Aa|Aa) -> (Aa) and (Aa)
	idx_1 = rand(1:ploidy(a))
	m1 = a.loci[idx_1, :]
	# 4. meiosis 2:   (Aa) -> (A) (a), (aa) -> (a) (a), (AA) -> (A) (A)
	l = length(a)
	g1 = m1[1:Int(l/2)]
	g2 = m1[Int(l/2)+1:Int(l)]
	g = [g1,g2]
	idx_2 = rand(1:ploidy(a))
	m2 = g[idx_2]
    return m2
end
# XXX: There seems to be some confusion on gamete formation. We should have
# recombination before gamete formation, and random segregation of homologs
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
	gametogenesis(d::AbstractDeme{A})

Gamete formation in a mixed ploidy population (1n to 4n).
"""
function gametogenesis(d::AbstractDeme{A}) where A
	new_agents = A[]
	for agent in d.agents
		push!(new_agents, gametogenesis(agent,d))
	end
    d(new_agents)
end

"""
	mateh(a::Agent, b::Agent)
Mating in a haploid population.
"""
function mateh(a::Agent, b::Agent)
    newloci = similar(a.loci)
    for i in 1:length(a)
		@inbounds newloci[i] = rand() < 0.5 ? a[i] : b[i]
    end
	Agent(newloci, a.d)
end	

"""
	random_matingh(d::AbstractDeme)
Random mating in a haploid deme.
"""
function random_matingh(d::AbstractDeme)
    newdeme = similar(d.agents)
    for i=1:length(d)
		@inbounds newdeme[i] = mateh(rand(d, 2)...)
	end
	d(newdeme)
end 

"""
	mate_p(a::Agent, b::Agent)
"""
function mate_p(a, b)
    ar = recombination(a)
    br = recombination(b)                                                                                                    
    Agent(loci=[ar.loci[1:1,:]; br.loci[1:1,:]], d=a.d)
end
	
"""
	mate_p(a::Agent, b::Agent, d::AbstractDeme)

Mating in a mixed ploidy deme. Assumes that different cytotypes can be
compatible with a decrease in viability (cfr. OffspringViability matrix) (i.e.
when #individuals with different ploidy hybridize, they generate have a
probability p to generate no viable offspring).  Selfing is allowed without
cost. This influences the dynamics by including hybrid offspring that can
compete for space (and affects the malthusian fitness/density dependence of
selection).
"""	
function mate_p(a::Agent, b::Agent, d::AbstractDeme)
    #gamete formation
    # XXX we should have recombinatoin *before* gamete formation ?!
    ag = gametogenesis(recombination(a), d)
    bg = gametogenesis(recombination(b), d)
    #combine gametes and assign viability
    via = viability(ag, bg, d) 
    # not sure if returning 0 is the best idea in terms of type stability
    # perhaps you should return a mock-agent (for instance an agent with length
    # zero genome) when there is no viable offspring
	# if d is defined at agent level, how do we cope with it here? 
    return rand() < via ? Agent(loci=[ag.loci ; bg.loci], d=a.d) : 0
end

"""
	random_mating(d::AbstractDeme{A}) where A

Random mating in a mixed ploidy deme.
"""
function random_mating(d::AbstractDeme{A}) where A
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
	number_of_offspring(d::AbstractDeme, a::Agent)
"""
function number_of_offspring(d::MixedPloidyDeme, a::Agent)
    logw = malthusian_fitness(d, a)
    rand(Poisson(exp(logw)))
end

"""
	number_of_offspring(d::IslandDeme, a::Agent)
"""
function number_of_offspring(d::IslandDeme, a::Agent)
    logw = directional_selection(d, a)
    rand(Poisson(exp(logw))) #does this kind of demographic stochasticity make sense with directional selection? 
	#with a poisson distributed number of offspring the variance increases with increasing fitness (exp(logw))
	#rand(Normal(logw,1))
end

"""
	mating_PnB(d::MixedPloidyDeme{A})

Mating in a mixed ploidy deme with unreduced gamete formation and partner
choice weighted by fitness where every individual has all it's offspring with the same partner(cfr. PnB).
"""
function mating_PnB(d::MixedPloidyDeme{A}) where A
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

"""
	mating_PnB_x(d::MixedPloidyDeme{A})

Mating in a mixed ploidy deme with unreduced gamete formation and partner
choice weighted by fitness.
"""
function mating_PnB_x(d::AbstractDeme{A}) where A
	new_agents = A[]
	fitnesses = exp.(malthusian_fitness(d))
	for i=1:length(d)
		B1 = d.agents[i]
		noff = number_of_offspring(d, B1)
        Bs = sample(d.agents, weights(fitnesses), noff) 
        offspring = filter(!ismock, map(B2->mate_p(B1, B2, d), Bs))
        new_agents = vcat(new_agents, offspring)
    end
    d(new_agents)
end

"""
	mating_PnB_xh(d::MixedPloidyDeme{A})

Mating in a haploid deme with selection.
"""
function mating_PnB_xh(d::AbstractDeme{A}) where A
	new_agents = A[]
	fitnesses = exp.(malthusian_fitness(d))
	for i=1:length(d)
		B1 = d.agents[i]
		noff = number_of_offspring(d, B1)
        Bs = sample(d.agents, weights(fitnesses), noff) 
        offspring = filter(!ismock, map(B2->mateh(B1, B2), Bs))
        new_agents = vcat(new_agents, offspring)
    end
    d(new_agents)
end

"""
	malthusian_fitness(d::AbstractDeme,a::Agent)
Return the Malthusian fitness (density dependence and stabilizing selection) of a single agent in a deme.
"""
function malthusian_fitness(d::MixedPloidyDeme,a::Agent)
    N = length(d)
    z = trait(a)
    return d.rm*(1-(N/d.K))-((z-d.θ)^2)/(2*d.Vs)
end 

"""
	malthusian_fitness(d::AbstractDeme)
Return the Malthusian fitness (density dependence and stabilizing selection) of each agent in a deme.
"""
function malthusian_fitness(d::MixedPloidyDeme)
    N = length(d)
    fitnesses = Float64[]
	for agent in d.agents
		z = trait(agent)
    	f = d.rm*(1-(N/d.K))-((z-d.θ)^2)/(2*d.Vs)
		push!(fitnesses, f)
	end
	fitnesses
end

"""
	directional_selection(d::AbstractDeme,a::Agent)
"""
function directional_selection(d::IslandDeme,a::Agent)
    N = length(d)
    z = trait(a)
    return d.β*(z-d.θ)
end 

"""
	directional_selection(d::IslandDeme)
"""
function directional_selection(d::IslandDeme)
    N = length(d)
    fitnesses = Float64[]
	for agent in d.agents
		z = trait(agent)
    	f = d.β*(z-d.θ)
		push!(fitnesses, f)
	end
	return fitnesses
end

"""
	mating_PnB(d::IslandDeme{A})

Mating in a mixed ploidy islanddeme with unreduced gamete formation directional selection.
"""
function mating_PnB(d::IslandDeme{A}) where A
	new_agents =  A[]
	fitnesses = exp.(directional_selection(d))
	for i=1:length(d)
 		B1 = d.agents[i]
		noff = number_of_offspring(d,B1)
		B2 = sample(d.agents, weights(fitnesses))
		m = mate_p(B1,B2,d)
		if m != 0
			for c in 1:noff 
				push!(new_agents, m)
			end
		end
	end
    d(new_agents)
end 

#Habitat level

"""
	linear_gradient(Dm,b,θ)

Initialize a vector of phenotypic optima for each deme of the habibat according
to a linear gradient.
- `b` : Steepness of linear gradient
- `θ` : Phenotypic optimum in the central deme of the habitat
- `Dm` : Number of demes 
"""
function linear_gradient(b,θ,Dm) 
	#KK = [i*b for i in 0:Dm-1]
	#KK = [θ for i in 0:Dm-1]
	a = -(Dm/2)*b + θ #where Dm is number of demes and b is gradient
	KK = [(i*b)+a for i in 1:Dm]
    return KK
end

"""
	initiate_habitat(d::MixedPloidyDeme, gradient, p, α, L, N)

Aim should be to initiate a population for nd_s demes on a linear gradient
(with slope b) and with optimal genetic variance where one half of the genes
are adapted, meaning their clines take the form and spacing as assumed for the
deterministic model under linkage equilibrium.
- `d` : Empty deme
- `gradient` : Linear gradient
- `p` : 
- `α` : Allelic effect size
- `L` : Number of loci
- `N` : Number of starting individuals in central deme
"""
function initiate_habitat(d::MixedPloidyDeme, gradient, p, α, L, N)
    #ag = randagent_p(p, α, L, [0., 1., 0., 0.], 0)
	#had to change some things back here:
	#it wasn't giving the expeted output somehow but filling all the demes with agents instead of just the middle one
	#also has to initialize every deme with different θ so not sure if it works with the object syntax
    hab = Habitat(demes=[MixedPloidyDeme(agents=randagent_p(p, α, L, [0., 1., 0., 0.], 0),K=d.K,θ=i,rm=d.rm,Vs=d.Vs,α=d.α,μ=d.μ,OV=d.OV,UG=d.UG) for i in gradient])
	for a in randagent_p(p, α, L, [0., 1., 0., 0.], N)
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

"""
	gaussian_dispersal(h::Habitat,σ)
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

#Mutation:

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
	evolving_haploiddeme(d::AbstractDeme, ngen)

Simulate a single random mating deme with mixed ploidy.
"""
function evolving_haploiddeme(d::AbstractDeme, ngen; heterozygosities_p = heterozygosities_p, allelefreqs_p = allelefreqs_p, 
	trait_mean = trait_mean, pf = ploidy_freq, fta = f_trait_agents)
	het = [heterozygosities_p(d)]
	af = [allelefreqs_p(d)]
	tm = [trait_mean(d)]
	p1 = [ploidy_freq(d)[1]]
	fta = [f_trait_agents(d)]
		
	for n=1:ngen
		d = random_matingh(d)
		#d = unreduced_gamete(d)
		push!(het, heterozygosities_p(d))
		push!(af, allelefreqs_p(d))
		push!(tm, trait_mean(d))
		push!(p1, ploidy_freq(d)[1])
		push!(fta, f_trait_agents(d))
	end
	(het=het, af=af, tm=tm, deme=d, p1=p1, ngen=ngen, fta=fta)
end

"""
	evolving_neutraldeme(d::AbstractDeme, ngen)

Simulate a single random mating deme with mixed ploidy.
"""
function evolving_neutraldeme(d::AbstractDeme, ngen; heterozygosities_p = heterozygosities_p, 
	allelefreqs_p = allelefreqs_p, trait_mean = trait_mean, pf = ploidy_freq, fta = f_trait_agents)
	het = [heterozygosities_p(d)]
	af = [allelefreqs_p(d)]
	tm = [trait_mean(d)]
	p2 = [ploidy_freq(d)[2]]
	p3 = [ploidy_freq(d)[3]]
	p4 = [ploidy_freq(d)[4]]
	fta = [f_trait_agents(d)]
	
	for n=1:ngen
		d = random_mating(d)
		#d = unreduced_gamete(d)
		push!(het, heterozygosities_p(d))
		push!(af, allelefreqs_p(d))
		push!(tm, trait_mean(d))
		push!(p2, ploidy_freq(d)[2])
		push!(p3, ploidy_freq(d)[3])
		push!(p4, ploidy_freq(d)[4])
		push!(fta, f_trait_agents(d))
	end
	(het=het, af=af, tm=tm, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, fta=fta)
end

"""
	evolving_selectiondeme(d::AbstractDeme, ngen)

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.
"""
function evolving_selectiondeme(d::MixedPloidyDeme, ngen; 
	heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = 
	allelefreqs_p, pf = ploidy_freq, fta = f_trait_agents)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[2]]
	p3 = [ploidy_freq(d)[3]]
	p4 = [ploidy_freq(d)[4]]
	fta = [f_trait_agents(d)]
	
	for n=1:ngen
		d = mating_PnB_x(d)
		d = mutate(d) 
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[2])
		push!(p3, ploidy_freq(d)[3])
		push!(p4, ploidy_freq(d)[4])
		push!(fta, f_trait_agents(d))
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, het=het,tm=tm, af=af, fta=fta) 
end

"""
	evolving_selectiondemeh(d::AbstractDeme, ngen)

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.
"""
function evolving_selectiondemeh(d::MixedPloidyDeme, ngen; 
	heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = 
	allelefreqs_p, pf = ploidy_freq, fta = f_trait_agents)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[2]]
	p3 = [ploidy_freq(d)[3]]
	p4 = [ploidy_freq(d)[4]]
	fta = [f_trait_agents(d)]
	
	for n=1:ngen
		d = mating_PnB_xh(d)
		d = mutate(d) 
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[2])
		push!(p3, ploidy_freq(d)[3])
		push!(p4, ploidy_freq(d)[4])
		push!(fta, f_trait_agents(d))
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, het=het,tm=tm, af=af, fta=fta) 
end
	
"""
	evolving_ugdeme(d::AbstractDeme, ngen)

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.
"""
function evolving_ugdeme(d::MixedPloidyDeme, ngen; pf = ploidy_freq)

	pop = [length(d)]
	p2 = [ploidy_freq(d)[2]]
	p3 = [ploidy_freq(d)[3]]
	p4 = [ploidy_freq(d)[4]]
	
	for n=1:ngen
		d = mating_PnB(d)
		d = mutate(d) 
		push!(pop, length(d))
		push!(p2, ploidy_freq(d)[2])
		push!(p3, ploidy_freq(d)[3])
		push!(p4, ploidy_freq(d)[4])
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen) 
end

"""
	evolving_islanddeme(d::AbstractDeme, ngen)

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.
"""
function evolving_islanddeme(d::IslandDeme, ngen; 
	heterozygosities_p=heterozygosities_p, fit=directional_selection, trait_mean = trait_mean, allelefreqs_p = 
	allelefreqs_p, pf = ploidy_freq, fta = f_trait_agents)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[2]]
	p3 = [ploidy_freq(d)[3]]
	p4 = [ploidy_freq(d)[4]]
	fta = [f_trait_agents(d)]
	
	for n=1:ngen
		d = mating_PnB(d)
		d = mutate(d) 
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[2])
		push!(p3, ploidy_freq(d)[3])
		push!(p4, ploidy_freq(d)[4])
		push!(fta, f_trait_agents(d))
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, het=het,tm=tm, af=af, fta=fta) 
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
"""
ploidy_freq(d::AbstractDeme) = counts(map(ploidy, d.agents), 1:4)

var_add(a::Agent,α) = ploidy(a)*α^2*sum(heterozygosities_p(a))

#Genotype -> phenotype maps

"""
    trait(a::Agent)
"""
trait(a::Agent) = sum(a)/ploidy(a)
	
"""
	trait_mean(d::AbstractDeme)
"""
function trait_mean(d::AbstractDeme)
    length(d.agents) == 0 && return 0.
    z = map(trait, d.agents)
    sum(z)/length(d)
end

"""
	f_trait_agents(d::AbstractDeme)
"""
f_trait_agents(d::AbstractDeme) = map(trait, d.agents)

"""
	loci_df(d::AbstractDeme)
"""
function loci_df(d::AbstractDeme)
	v = [agent.loci for agent in d.agents]
	h = reduce(vcat,v)
	DataFrame(h)
end

#dominance
#pollen <-> seed migration
#assortative mating