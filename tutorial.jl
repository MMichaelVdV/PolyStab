### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 1e540d75-0a9d-4bb2-85aa-02ea78f08766
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, DataFrames, StatsPlots

# ╔═╡ 42a769f0-1bf8-11ec-2a76-e7cea7f2e77a
import PolyStab:randagent

# ╔═╡ 9fe6dd08-89c2-49c7-99ed-bf27e4a90242
md"##### Generating agents"

# ╔═╡ 111bd54b-d480-4c5e-9a92-aa6f38a2d074
md""" Generating a single agent:"""

# ╔═╡ ee0a9506-823f-449d-8bb9-7f864c8e73ae
"""
    randagent(p, α, n, k; d)

-n: number of loci
-k: ploidy l
-α: allelic effect size
-p: starting allele frequency
"""
a = randagent(0.5, 0.1, 100, 2, d=1.)

# ╔═╡ 54ea2331-7430-48c6-b0b9-f82c4b80f3ff
a.loci

# ╔═╡ 1c163c4d-d8dd-49a4-a07a-372abf74ef21
import PolyStab:trait

# ╔═╡ 3f95bac8-5c98-472f-b53d-f5f2ebe2e113
trait(a)

# ╔═╡ 2e17fd8e-f360-421e-ba35-8921c9807c47
md""" It is also possible to generate an array of agents:"""

# ╔═╡ 4d4f91c4-116e-45bd-aca5-a4e728026159
import PolyStab:randagent_p

# ╔═╡ de22de55-116d-41fd-8865-a9a0ae5be0c9
"""
    randagent_P(p, α, n, k, N; d)

-n: number of loci
-k: ploidy l
-α: allelic effect size
-p: starting allele frequency
-N: number of individuals
"""
b = randagent_p(0.5, 0.1, 100, [0.,1.,0.,0.], 10, d=1.)

# ╔═╡ 7c3af315-31b4-40bb-af86-8a105f427b7c
b

# ╔═╡ 5e0aaf09-4453-4aea-b0f9-5f14b16200f5
md"""#### Generating a single mixed-ploidy deme"""

# ╔═╡ 870aab66-3e7a-4b47-917d-190a72f8e2de
import PolyStab:MixedPloidyDeme

# ╔═╡ 5d5393dd-bebe-4275-992b-7477a5c7834b
"""
    MixedPloidyDeme{A,T}

Single mixed-ploidy deme, most of the 'population
genetic' environment is implemented at this level (drift, selection,
mutation). OV: Viability matrix, a symmetric matrix that contains the viability
of offspring for each possible combination of gametes. Ug: Unreduced gamete
formation rate matrix, a matrix that contains the probability of unreduced gamete formation for each level of ploidy in the population.

- `K` : Carrying capacity
- `θ` : Environmental optimum
- `rm`: Growth rate per capita
- `Vs`: Variance of stabilizing selection
- `α` : Allelic effect size
- `μ` : Mutation rate
- `OV` : Offspring viability
- `UG` : Unreduced gamete probabilities
"""
d_p1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [1., 0., 0., 0.],10), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 0.1)

# ╔═╡ 6b866b78-80d7-43b1-a4c5-33f6c7908a03
d_p1

# ╔═╡ 1dc54ec2-f4d3-4237-8fb2-b331e0bf887f
d_p1.agents

# ╔═╡ 6e07be7d-8af6-4bf1-9591-f42ab3ea9db0
trait.(d_p1.agents)

# ╔═╡ a24f94e4-c012-4c82-ba4c-93484d0a7263
md"""#### Simulating a single mixed-ploidy deme"""

# ╔═╡ e2ced766-0032-4cb2-b73c-fd1a086e0950
import PolyStab:evolving_selectiondeme

# ╔═╡ 7ac2f20c-db6a-4927-a76b-d7ded656adbc
"""
	evolving_selectiondeme(d::AbstractDeme, ngen)

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.

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


# ╔═╡ bb6aa892-6ef1-480a-8a2e-be1c131c0d5b
evolving_selectiondeme(d_p1, 500)

# ╔═╡ f91908fb-8202-4087-9ee0-ed8207a41378
md""" ###### Components of the mating function"""

# ╔═╡ a50b4da8-39d5-4977-a1de-a8a172c0e6ac
"""
	mating_PnB_x(d::MixedPloidyDeme{A})

Mating in a mixed ploidy deme with unreduced gamete formation and partner
choice weighted by malthusian fitness.

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

# ╔═╡ d6eeb3b0-fbce-4ed8-b32b-0adbae74019c
"""
	malthusian_fitness(d::AbstractDeme)
Return the Malthusian fitness (density dependence and stabilizing selection) of each agent in a deme.

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

# ╔═╡ 3db37b1c-693b-4c34-9e18-b9126b0fe235
"""
	number_of_offspring(d::AbstractDeme, a::Agent)

function number_of_offspring(d::MixedPloidyDeme, a::Agent)
    logw = malthusian_fitness(d, a)
    rand(Poisson(exp(logw)))
end
"""

# ╔═╡ 776a2d51-ca39-4368-9431-5a3ad8afe14d
"""
	mate_p(a::Agent, b::Agent, d::AbstractDeme)

Mating in a mixed ploidy deme. Assumes that different cytotypes can be
compatible with a decrease in viability (cfr. OffspringViability matrix) (i.e.
when #individuals with different ploidy hybridize, they generate have a
probability p to generate no viable offspring).  Selfing is allowed without
cost. This influences the dynamics by including hybrid offspring that can
compete for space (and affects the malthusian fitness/density dependence of
selection).
	
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

# ╔═╡ 82e7cf93-ba18-483a-9880-4dd396416510
"""
	recombination(a::Agent)

Free recombination between loci in a mixed ploidy population. 

function recombination(a::Agent)
    genome = similar(a.loci)
    for i in 1:nloci(a)
        genome[:,i] = shuffle(a.loci[:,i])
	end
    Agent(loci=genome, d=a.d) 
end


	gametogenesis(a::Agent, d::AbstractDeme)

Gamete formation in a mixed ploidy population (1n to 4n).

function gametogenesis(a::Agent, d::AbstractDeme)
	num = sample(1:4, weights(d.UG[ploidy(a),:]))
    idx = sample(1:ploidy(a), num, replace=false)
    return Agent(a.loci[idx, :], a.d)
end	


	viability(a::Agent, b::Agent, d::deme)
After fusing of gametes, calculate the viability of resulting individual. For example the viability of
a triploid individual is by default set to 0, diploid and tetraploid individuals set to 1.

function viability(a::Agent, b::Agent, d::AbstractDeme)
	return d.OV[ploidy(a), ploidy(b)]
end
"""


# ╔═╡ Cell order:
# ╠═1e540d75-0a9d-4bb2-85aa-02ea78f08766
# ╠═42a769f0-1bf8-11ec-2a76-e7cea7f2e77a
# ╟─9fe6dd08-89c2-49c7-99ed-bf27e4a90242
# ╟─111bd54b-d480-4c5e-9a92-aa6f38a2d074
# ╠═ee0a9506-823f-449d-8bb9-7f864c8e73ae
# ╠═54ea2331-7430-48c6-b0b9-f82c4b80f3ff
# ╠═1c163c4d-d8dd-49a4-a07a-372abf74ef21
# ╠═3f95bac8-5c98-472f-b53d-f5f2ebe2e113
# ╟─2e17fd8e-f360-421e-ba35-8921c9807c47
# ╠═4d4f91c4-116e-45bd-aca5-a4e728026159
# ╠═de22de55-116d-41fd-8865-a9a0ae5be0c9
# ╠═7c3af315-31b4-40bb-af86-8a105f427b7c
# ╟─5e0aaf09-4453-4aea-b0f9-5f14b16200f5
# ╠═870aab66-3e7a-4b47-917d-190a72f8e2de
# ╠═5d5393dd-bebe-4275-992b-7477a5c7834b
# ╠═6b866b78-80d7-43b1-a4c5-33f6c7908a03
# ╠═1dc54ec2-f4d3-4237-8fb2-b331e0bf887f
# ╠═6e07be7d-8af6-4bf1-9591-f42ab3ea9db0
# ╟─a24f94e4-c012-4c82-ba4c-93484d0a7263
# ╠═e2ced766-0032-4cb2-b73c-fd1a086e0950
# ╠═7ac2f20c-db6a-4927-a76b-d7ded656adbc
# ╠═bb6aa892-6ef1-480a-8a2e-be1c131c0d5b
# ╟─f91908fb-8202-4087-9ee0-ed8207a41378
# ╠═a50b4da8-39d5-4977-a1de-a8a172c0e6ac
# ╠═d6eeb3b0-fbce-4ed8-b32b-0adbae74019c
# ╠═3db37b1c-693b-4c34-9e18-b9126b0fe235
# ╠═776a2d51-ca39-4368-9431-5a3ad8afe14d
# ╠═82e7cf93-ba18-483a-9880-4dd396416510
