### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ c0283090-882c-11eb-0f37-c546cb7cb8b9
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes

# ╔═╡ cf15a560-8822-11eb-0789-ad5bb86f481b
md""" Some testing to transfer from array to matrix"""

# ╔═╡ 38cd6f40-882f-11eb-3416-1fdfde76f20d
begin

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
@with_kw struct MixedPloidyDeme{A,T} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 100
    θ ::T     = 12.5
    rm::T     = 1.06
    Vs::T     = 1/2
    u ::T     = 0.01
    μ ::T     = 1e-6
	OV::Matrix{T} = [1. 0.1 0. 0. ; 0.1 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.]
	UG::Matrix{T} = [0.9 0.1 0. 0. ; 0. 0. 0. 0. ; 0. 0.9 0. 0.1]
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
	
end

# ╔═╡ 1f6f2ebe-8830-11eb-159f-ef99983ae560
d_p = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 0.5, 50, [2], 50, d=2.),randagent_p(0.5, 0.5, 50, [4], 0, d=4.)), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0.9 0.1 0. 0. ; 0. 0. 0. 0. ; 0. 0.9 0. 0.1] )

# ╔═╡ ab91dec0-8835-11eb-2594-bb2fed4f14e7
a = randagent_p(0.5, 0.5, 50, 4, d=4.)

# ╔═╡ 72920d0e-8908-11eb-1fc7-7f2ec0184880
a1 = randagent_p(0.5, 0.5, 50, 1, d=1.)

# ╔═╡ 72a93e90-8908-11eb-37f8-cd6ad67826a4
a2 = randagent_p(0.5, 0.5, 50, 2, d=2.)

# ╔═╡ ebcf3950-8836-11eb-1092-9fe734282054
d_p.agents

# ╔═╡ b86086ae-8835-11eb-225b-bf9c29114e21
size(a.loci)

# ╔═╡ 58da38c0-8836-11eb-18c7-e3087770b4fb
length(a)

# ╔═╡ c127ec50-8837-11eb-21d1-0b4c07df9464
sum(a)

# ╔═╡ d21432d0-8837-11eb-11bc-49a4315c59ff
sum(a.loci)

# ╔═╡ eb74ab50-8838-11eb-25e5-010384ccca19
d_p.UG[3,:]

# ╔═╡ c636e860-883a-11eb-27e7-c9e3f70fbe75
num = sample([1.,2.,3.,4.],weights(d_p.UG[ploidy(a)-1,:]))

# ╔═╡ 324f5502-883b-11eb-15de-4fd9dbdbe959
loci = (Float64[], num)

# ╔═╡ cab6126e-883b-11eb-2fc4-65f85a389ffd
a.loci[1,:]

# ╔═╡ 3c6ace60-883c-11eb-215f-f95e755463d9
a.loci[shuffle(1:end), :]

# ╔═╡ fcb619e0-883c-11eb-3bae-47a43a41f4d3
begin
	l = zeros(4, length(a))
	b = deepcopy(a.loci)
	i = b[1,:]
	b = b[1:end .!= 2,:]
end

# ╔═╡ 3ee1cc60-883d-11eb-041c-3b49db5465cf
l

# ╔═╡ aa77f820-88d6-11eb-1f29-959c712ce532
[a.loci ; a.loci][1,1]

# ╔═╡ 2a2d3dce-883e-11eb-2c87-4343cac53699
similar(a.loci[1,:])

# ╔═╡ 34fecc50-8835-11eb-0e1f-4b2e13ff740e
begin
	
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
function recombine_poly(d::AbstractDeme{A}) where A
	new_agents =  A[]
	for agent in d.agents
		num = ploidy(agent)
		loci = zeros(num, length(a))
		newlocus = similar(agent.loci[1,:])
		
		for l in 1:num
			for j in 1:length(a) #this loops over the different loci for each chrosome
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
	num = sample([1.,2.,3.,4.],weights(d.UG[ploidy(a)-1,:]))
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

number_of_offspring(d::AbstractDeme,a::Agent) = rand(Poisson(exp(malthusian_fitness(d::AbstractDeme,a::Agent))))

"""
"""
trait(a::Agent) = sum(a)/a.d
	
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
	heterozygosities_p(d::AbstractDeme)
"""
function heterozygosities_p(d::AbstractDeme)
	freqs=allelefreqs_p(d)
    map(p->p*(1-p), freqs)
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

end

# ╔═╡ 4fd07fa0-8908-11eb-11aa-c5f5c06cec55
viability(a1, a2, d_p) 

# ╔═╡ c6c83340-8837-11eb-048c-f91ba08f9826
trait(a)

# ╔═╡ 0e650dc0-8835-11eb-122e-f1b02f2471a0
sim_ploidyvar1 = evolving_deme_ploidyvar(d_p,50)

# ╔═╡ 09055fc0-88d9-11eb-27ea-032e5733ac17
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

# ╔═╡ c3953840-8836-11eb-2b8b-516d9074e45d
ploidy_freq(d_p)

# ╔═╡ da90493e-88d6-11eb-1873-8127ca6a1650
begin
	#gamete formation
	
	ag = recombine_poly(unreduced_gamete(a,d_p))
	bg = recombine_poly(unreduced_gamete(a,d_p))
	#combine gametes and assign viability
	via = viability(ag, bg, d_p) 
	rand() < via ? Agent([ag.loci ; bg.loci], 1. *(ploidy(ag) + ploidy(bg))) : 0
	
end

# ╔═╡ d15eb34e-88d8-11eb-145a-0f40e5907c42
begin
	recombine_poly(d_p)
end

# ╔═╡ f6e3deb0-8906-11eb-273d-b127626fa404
mate_p(a,a,d_p)

# ╔═╡ Cell order:
# ╟─cf15a560-8822-11eb-0789-ad5bb86f481b
# ╠═c0283090-882c-11eb-0f37-c546cb7cb8b9
# ╠═38cd6f40-882f-11eb-3416-1fdfde76f20d
# ╠═1f6f2ebe-8830-11eb-159f-ef99983ae560
# ╠═ab91dec0-8835-11eb-2594-bb2fed4f14e7
# ╠═72920d0e-8908-11eb-1fc7-7f2ec0184880
# ╠═72a93e90-8908-11eb-37f8-cd6ad67826a4
# ╠═4fd07fa0-8908-11eb-11aa-c5f5c06cec55
# ╠═ebcf3950-8836-11eb-1092-9fe734282054
# ╠═b86086ae-8835-11eb-225b-bf9c29114e21
# ╠═58da38c0-8836-11eb-18c7-e3087770b4fb
# ╠═c127ec50-8837-11eb-21d1-0b4c07df9464
# ╠═c6c83340-8837-11eb-048c-f91ba08f9826
# ╠═d21432d0-8837-11eb-11bc-49a4315c59ff
# ╠═0e650dc0-8835-11eb-122e-f1b02f2471a0
# ╠═09055fc0-88d9-11eb-27ea-032e5733ac17
# ╠═c3953840-8836-11eb-2b8b-516d9074e45d
# ╠═eb74ab50-8838-11eb-25e5-010384ccca19
# ╠═c636e860-883a-11eb-27e7-c9e3f70fbe75
# ╠═324f5502-883b-11eb-15de-4fd9dbdbe959
# ╠═cab6126e-883b-11eb-2fc4-65f85a389ffd
# ╠═3c6ace60-883c-11eb-215f-f95e755463d9
# ╠═fcb619e0-883c-11eb-3bae-47a43a41f4d3
# ╠═3ee1cc60-883d-11eb-041c-3b49db5465cf
# ╠═aa77f820-88d6-11eb-1f29-959c712ce532
# ╠═2a2d3dce-883e-11eb-2c87-4343cac53699
# ╠═da90493e-88d6-11eb-1873-8127ca6a1650
# ╠═d15eb34e-88d8-11eb-145a-0f40e5907c42
# ╠═f6e3deb0-8906-11eb-273d-b127626fa404
# ╠═34fecc50-8835-11eb-0e1f-4b2e13ff740e
