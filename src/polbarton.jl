# Think of how we will represent polyploids. We will need to hold different
# chromosomes, not sure if a vector of vectors or a matrix is best. Also we
# will need to be able to scale allelic effects, using a parameter `d`. The
# latter only comes in when we compute the phenotype.
# Should the mutation rate and unreduced gamete formation rate be individual
# level parameters? Probably they should...
"""
    Agent{T}

Should be an agent of arbitrary (allowed) ploidy level.
"""
@with_kw struct Agent{T}
    loci::Vector{T}
    d::T = 1. # allelic effect scaler for different ploidy levels
end

Base.length(a::Agent) = length(a.loci)

# generalize to [i,j] indexing for different ploidy levels
Base.getindex(a::Agent, i) = a.loci[i]  

# generalize to different ploidy levels
randagent(p, α, n; d=1.) = Agent(rand(Bernoulli(p), n) * α, d)
randagent(p, α, n, N; d=1.) = [randagent(p, α, n, d=d) for i=1:N]

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


# Note, I think it makes sense to have an `AbstractDeme` type, because we can
# easily imagine that when we want to implement more complex models (say with
# different fertility across ploidy levels or assortative mating or whatever)
# that we implement a new `Deme` type, but call the same functions like
# `random_mating` or whatever. So we can try to implement the main simulation
# of a habitat for the `AbstractDeme` type so that it works for any concrete 
# deme type that implements the necesary functions (such as `random_mating`)
abstract type AbstractDeme{A} end

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

# construct a new deme from an old one, this uses the 'type as function' syntax
(d::MixedPloidyDeme)(agents) = MixedPloidyDeme(agents, d.K, d.θ, d.rm, d.Vs, d.u, d.μ)

@with_kw struct SimpleDeme{A} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 50
end

(d::SimpleDeme)(agents) = SimpleDeme(agents, d.K)

Base.length(d::AbstractDeme) = length(d.agents)
Base.getindex(d::AbstractDeme, i) = d.agents[i]
Base.rand(rng::AbstractRNG, d::AbstractDeme, n) = rand(rng, d.agents, n)
Base.rand(d::AbstractDeme, n::Int) = rand(default_rng(), d.agents, n)


# implemented for abstract deme type
function random_mating(d::AbstractDeme)
    newdeme = similar(d.agents)
    for i=1:length(d)
		@inbounds newdeme[i] = mate(rand(d, 2)...)
	end
	d(newdeme)
end 

# won't work for mixed-ploidy
function allelefreqs(d::AbstractDeme)
    # assuming 0 is absence of allele
    f = j->mapreduce(i->d[i][j] != 0., +, 1:length(d))/length(d)
	frequencies = map(f, 1:length(d[1]))
end

# assumes haploid
function heterozygosities(d::AbstractDeme, freqs=allelefreq(d)) 
    map(p->p*(1-p), freqs)
end

"""
    Habitat{D}

A 1-dimensional habitat, i.e. an array of connected demes. This implements
the migration aspects of the population genetic environment.
"""
struct Habitat{D,T}
    demes::Vector{D}
    σ ::T
    b ::T
    Dm::T
end
