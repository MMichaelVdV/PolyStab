### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ cc1dcf70-7c6c-11eb-3664-b7e70e49723a
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes

# ╔═╡ e419ec50-87fa-11eb-230a-8f6396f0e459
using Turing, MCMCChains, StatsPlots, DataFrames

# ╔═╡ 0144f6d0-87fb-11eb-0fa7-7149bfb316db
using StatsFuns: logistic

# ╔═╡ bd2ab750-7c60-11eb-3cc5-6d57dd448dfd
md"""##### Minority Cytotype Exclusion in Local Plant Populations (Levin, 1995)"""

# ╔═╡ 98180cf0-7c6b-11eb-0da9-c55b3cf2a066


# ╔═╡ 2dadbbf0-7c5f-11eb-3d09-1fcf2b82ead9
md"""##### Establishment of a tetraploid cytotype in a diploid population: Effect of relative fitness of the cytotypes (Felber, 1991)"""

# ╔═╡ 472102d2-84a8-11eb-2d6d-eb101b7a707b


# ╔═╡ 7df71be0-7c70-11eb-3b92-91457ecf38fd
md"""### Case 1 (Felber, 1991)"""

# ╔═╡ c8734230-7c65-11eb-30bd-4314f44ee3ba
md"""  
- Fertility of diploids and tetraploids is the same
- Viability of diploids and tetraploids is the same
- Triploids are not viable
- (The system is at equilibrium)

The only parameter in this model is u: unreduced gamete formation for diploids.

Starting from a population of only diploid individuals, according to this model, tetraploids will exclude diploids if unreduced gamete formation for is higher than 0.1716.

"""

# ╔═╡ 1f6a4e72-7d48-11eb-0933-eb25e7a5b568
md"""I used grid search on an interval for from 0 to 0.5 with 200 steps. For each value of u I simulated 50 generations and checked the most frequent ploidy level at the final generation"""

# ╔═╡ d75aee00-7c84-11eb-2d5a-b3ca3c4bd625
function bin(data)
	up2 = []
	up4 = []
	for i in data
		
		if i[2] == 2
			push!(up2,i[1])
		else
			push!(up4,i[1])
		end
		
	end
	up2,up4
end

# ╔═╡ a478224e-7c70-11eb-133b-03e027bf81c8
md""" ### Case 2 (Felber, 1991)"""

# ╔═╡ 57d84bf0-7c6b-11eb-0e68-f99c717f46e4
md"""  
See case 1 but now the parameter in this model is u(=v): unreduced gamete formation for both diploids and tetraploids.

Starting from a population of only diploid individuals, according to this model, tetraploids will exclude diploids if unreduced gamete formation for is higher than 0.2.

"""

# ╔═╡ 1abbabe0-87fb-11eb-271e-53e02f1ec547
md""" ### Cytogenic load"""

# ╔═╡ 86be28f0-84ad-11eb-2ccf-6f60ad079d75
begin
#credits to arzwa
#This is the simple Felber model, expressed in terms of marginal/relative fitnesses
w = [0., 1., 0., 1.]
# wij = proportion of offspring of ploidy j for a given individual of ploidy i
# marginal/relative fitness of ploidy i (proportion of offspring in the next generation coming from a parent of ploidy i, regardless of offspring ploidy level) is
# wi* = ∑ⱼwij. The mean fitness is then ∑ᵢwi* pᵢ. The proportion of offspring
# of ploidy level j in the next generation is (∑ᵢwij pᵢ) / w̄. This seems like a
# mathematically quite elegant formulation of the model, and one that may
# generalize relatively easily.
w22(w, u, p) = w[2] * (1-u)^2 * p
w24(w, u, p) = w[4] * (u^2 * p + u*(1-u) * (1-p))
w44(w, u, p) = w[4] * ((1-u)^2 * (1-p))+ u*(1-u) * p #u*(1-u) * p: where does this come from ?
w̄(w, u, p) = w22(w, u, p) * p  + w44(w, u, p) * (1-p) + w24(w, u, p) * p     
function w̄_(w, u, p) 
    w2 = w22(w, u, p) * p  
    w4 = w44(w, u, p) * (1-p) + w24(w, u, p) * p 
    return w2, w4
end
function evolve(w, u, p, n)
    map(1:n) do x
        w2, w4 = w̄_(w, u, p)
        p = w2 / (w2 + w4)
        p, w2 + w4
    end
end
# In this model (no triploids, tetraploids and diploids equally fit a priori,
# shared unreduced gamete formation rate, no higher polyploids) tetraploids
# take over irrespective of the initial frequencies as soon as u > 0.2. Note
# though that values of u > 0.2 lead to a very serious load in the population,
# as at least a fraction (1-0.8^2) of inviable offspring will be generated
# each generation!

p = plot(grid=false, legend=false, xlabel="\$t\$", ylabel="\$\\bar{w}\$")
cols = cgrad(:magma, 0.05:0.01:0.45)
for u=0.05:0.01:0.45
    mw = last.(evolve(w, u, 1., 30))
    plot!(p, mw, color=cols[u], linewidth=2, alpha=0.5)
end
plot(p, size=(750,750))
#savefig("felber-wbar.pdf")
# Interestingly, mean fitness decreases in this system. Of course it is well
# known that in systems with frequency dependence mean fitness does not
# necessarily increase, but somehow I never checked this for the Felber model.
end

# ╔═╡ 7f32ce70-85ab-11eb-1c30-1ffd857a2545
evolve([0., 1., 0., 1.], 0.2, 1., 20)

# ╔═╡ f641bd10-85aa-11eb-03fc-0de54af93d35
begin
pf = plot(grid=false, legend=true, xlabel="\$t\$", ylabel="\$\\bar{w}\$")
#cols = cgrad(:magma, 0.05:0.01:0.45)
#for u=0.05:0.01:0.45
	u = 0.25
	s = (evolve(w, u, 1., 30))
    mw = last.(s)
	pfs = first.(s)
	#f = 
    plot!(pf, mw, color=cols[u], linewidth=2, alpha=0.5, label="mean fitness")
	plot!(pf, pfs, color="blue", linewidth=2, alpha=0.5, label="diploid freq")
	plot!(pf, 1 .-pfs, color="red", linewidth=2, alpha=0.5, label="tetrploid freq")
#end
plot(pf, size=(750,750))
#savefig("felber-fubar.pdf")
# Interestingly, mean fitness decreases in this system. Of course it is well
# known that in systems with frequency dependence mean fitness does not
# necessarily increase, but somehow I never checked this for the Felber model.
end

# ╔═╡ 8a935580-85ad-11eb-1c01-fbfe6528c0e2
(evolve(w, u, 1., 30))

# ╔═╡ 495309e0-85ac-11eb-3231-b3dce7ae61c5
first.(evolve([0., 1., 0., 1.], 0.5, 1., 30))

# ╔═╡ aab09870-84b1-11eb-11c6-ffa1224114d8
begin
#Here I try to plot the cytogenic load for different levels of u to check if there is indeed a maximimum for a mixed population
cytoload_2n(u) = 2*(1-u)*u
cytoload_4n(u) = 2*(1-u)*u+u^2
#cytoload_total(u,f) = (2*(1-u)*u)*(1-f) + (2*(1-u)*u+u^2)*(f)
#as total frequency of unviable gametes
cytoload_total(u,f) = (1-u)*f*(1-u)*(1-f)+(1-u)*f*u*(1-f)+(1-u)*f*u*f+(1-u)*(1-f)*u*(1-f)+(u*(1-f))^2+u*f*u*(1-f)
	
f=[i for i in 0.0:0.01:1.]

pl = plot(grid=false, legend=false, xlabel="\$f(4N)\$", ylabel="\$cytoload\$")
#cols = cgrad(:magma, 0.05:0.01:0.45)
for u=0.05:0.01:0.45
	load = cytoload_total.(u,f)
    plot!(f, load, color=cols[u], linewidth=2, alpha=0.5)
end
plot(pl, size=(750,750))
	
end

# ╔═╡ 302363ee-7c6c-11eb-0d7b-effeadcdeb60
md"""#### Functions"""

# ╔═╡ 41a5cff0-7c6c-11eb-2b03-d7337d430df4
begin

"""
    Agent{T}
Should be an agent of arbitrary (allowed) ploidy level.
"""
@with_kw struct Agent{T,N}
    loci::Array{Array{T,1},N}
    d::T = 1. # allelic effect scaler for different ploidy levels
end

abstract type AbstractDeme{A} 
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

@with_kw struct SimpleDeme{A} <: AbstractDeme{A}
    agents::Vector{A}
    K ::Int64 = 15
end
	
#Viability matrix, should probably become integrated within deme struct
	
struct OffspringViability{T,N}
    viability::Array{Array{T,1},N}
end

#Unreduced gamete probability matrix, should probably become integrated within deme struct
	
struct UnreducedGamete{T,N} 
    prob::Array{Array{T,1},N}
end
	
# construct a new deme from an old one, this uses the 'type as function' syntax
(d::MixedPloidyDeme)(agents) = MixedPloidyDeme(agents, d.K, d.θ, d.rm, d.Vs, d.u, d.μ)
(d::SimpleDeme)(agents) = SimpleDeme(agents, d.K)




	
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
Base.length(d::AbstractDeme) = length(d.agents)
Base.getindex(d::AbstractDeme, i) = d.agents[i]
Base.rand(rng::AbstractRNG, d::AbstractDeme, n) = rand(rng, d.agents, n)
Base.rand(d::AbstractDeme, n::Int) = rand(d.agents, n)
Base.push!(d::AbstractDeme, a::Agent) = push!(d.agents, a)

function mate(a::Agent, b::Agent)
    newloci = similar(a.loci)
    for i in 1:length(a)
		@inbounds newloci[i] = rand() < 0.5 ? a[i] : b[i]
    end
	Agent(newloci, a.d)
end	

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
	
function allelefreqs(d::AbstractDeme)
    # assuming 0 is absence of allele
    f = j->mapreduce(i->d[i][j] != 0., +, 1:length(d))/length(d)
	frequencies = map(f, 1:length(d[1]))
end

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

function heterozygosities(d::AbstractDeme) 
	freqs=allelefreq(d)
    map(p->p*(1-p), freqs)
end

function heterozygosities_p(d::AbstractDeme)
	freqs=allelefreqs_p(d)
    map(p->p*(1-p), freqs)
end

expected_heterozygosity(H₀, t, N) = ((1.0-1.0/N)^t)*H₀

Base.sum(a::Agent) = sum(sum(a.loci))
	
function trait(a::Agent)
	sum(a)/a.d
end

function trait_mean(d::AbstractDeme)
	z = Float64[]
	for agent in d.agents
		push!(z,trait(agent))
	end
	sum(z)/length(d)
end

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

ploidy(d::AbstractDeme) = [length(a.loci) for a in d.agents]

function unreduced_gamete(d::AbstractDeme{A}) where A
	new_agents =  A[]
	for agent in d.agents
		push!(new_agents, unreduced_gamete(d,agent))
	end
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end	

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

function malthusian_fitness(d::AbstractDeme,a::Agent)
    N = length(d)
    z = trait(a)
    return d.rm*(1-(N/d.K))-((z-d.θ)^2)/(2*d.Vs)
end 

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

function viability(a::Agent, b::Agent, OV::OffspringViability)
	return OV.viability[ploidy(a)][ploidy(b)]
end

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

function unreduced_gamete(d::AbstractDeme{A}, UG::UnreducedGamete) where A
	new_agents =  A[]
	for agent in d.agents
		push!(new_agents, unreduced_gamete(agent,UG))
	end
	MixedPloidyDeme(agents=new_agents,k=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end		
	
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
	MixedPloidyDeme(agents=new_agents)
end	

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
	
	
end

# ╔═╡ 8f919600-7c6b-11eb-33b7-bd5a762fb3eb
#This is a mixed ploidy deme startin with only diploids (p,α,number of alleles,ploidy level,number of individuals,phentonype scaler)
d_p = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 0.5, 45, [2], 175, d=2.),randagent_p(0.5, 0.5, 45, [4], 0, d=4.)), u=0.01, K=200)

# ╔═╡ e0789de0-7c6e-11eb-0900-3387c904ea0b
sim_popvar = evolving_deme_popvar(d_p,150)

# ╔═╡ e5916552-7c6e-11eb-1868-a7b856e0ad1f
begin
	pf1_var = sim_popvar.p2
	pf2_var = sim_popvar.p4
	plot(pf1_var, grid=false, color=:blue, label="diploid", title="u=0.01")
	plot!(pf2_var, grid=false, color=:red, label="tetraploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ 1d6f4e60-7c6f-11eb-0be9-9341661ea1d2
#Allocates probabalities for unreduced gamete formation for 2n, 3n, 4n individuals respectively, for example [0.9,0.1,0.,0.] -> 2n ind has probability of 0.9 to produce 1n gametes, 0.1 for 2n (unreduced) gamete,...this doesn't take into account aneuploid gametes yet.
UG_1 = UnreducedGamete([[0.9,0.1,0.,0.],[0.,0.,0.,0.],[0.0,1.,0.0,0.]])

# ╔═╡ 20380b50-7c6f-11eb-14c9-fb3fc36e318b
#On the rows/columns (symmetric matrix) are the ploidy levels of the gametes to be combined (1n to 4n).The upper bound is now at 4n, all ploidy levels above are not viable (ie. combination of 2n and 3n gamete has 0 viability).
OV_1 = OffspringViability([[1.,0.,0.1,0.],[0.,1.,0.,0.],[0.1,0.,0.,0.],[0.,0.,0.,0.]])

# ╔═╡ 8caaa5de-7c6f-11eb-387a-05cc2add5483
sim_ploidyvar1 = evolving_deme_ploidyvar(d_p,500,UG_1,OV_1)

# ╔═╡ 13d994e0-7c70-11eb-1617-c37054ccc63b
begin
	pf2_p1 = sim_ploidyvar1.p2
	pf3_p1 = sim_ploidyvar1.p3
	pf4_p1 = sim_ploidyvar1.p4
	plot(pf2_p1, grid=false, color=:blue, label="diploid",title="u=0.01")
	#plot!(pf3_p1, grid=false, color=:green, label="triploid")
	plot!(pf4_p1, grid=false, color=:red, label="tetraploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ 4cde5ba0-7c79-11eb-23e1-d1ac938c0e6b
function random_search(t)
	OV = OffspringViability([[1.,0.,0.1,0.],[0.,1.,0.,0.],[0.1,0.,0.,0.],[0.,0.,0.,0.]])
	ploidy = []
	param = []
	for i in 1:t
		u = rand(Uniform(0., 0.5))
		UG = UnreducedGamete([[1-u,u,0.,0.],[0.,0.,0.,0.],[0.0,1.,0.0,0.]])
		sim_ploidyvar = evolving_deme_ploidyvar(d_p,100,UG,OV)
		if sim_ploidyvar.p2[end] >= sim_ploidyvar.p4[end]
			push!(ploidy,2)
		else
			push!(ploidy,4)
		end
		push!(param, u)
	end
	ploidy, param
end
		

# ╔═╡ e315b920-7c7d-11eb-0f34-5782accf7286
function grid_search(t)
	OV = OffspringViability([[1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
	ploidy = []
	param = []
	pop_size = []
	for u in range(0, stop=0.5, length=t)
		UG = UnreducedGamete([[1-u,u,0.,0.],[0.,0.,0.,0.],[0.0,1.,0.0,0.]])
		sim_ploidyvar = evolving_deme_ploidyvar(d_p,50,UG,OV)
		if sim_ploidyvar.p2[end] >= sim_ploidyvar.p4[end]
			push!(ploidy,2)
		else
			push!(ploidy,4)
		end
		push!(param, u)
		pop = (sim_ploidyvar.p2[end] + sim_ploidyvar.p4[end])
		push!(pop_size, pop)
	end
	ploidy, param, pop_size
end

# ╔═╡ 39c92860-7c7e-11eb-1efe-13bf99e99d1e
stats_2 = grid_search(200)

# ╔═╡ 5c9c61c0-87f2-11eb-107c-eb6817ebc523
begin
	M = [stats_2[1] stats_2[2]];
	data = DataFrame(M, ["Ploidy","u"]);
end

# ╔═╡ ec061ade-8801-11eb-2449-f3788b4ffa9b
begin
using CSV
CSV.write("./model_1.csv", data)
end

# ╔═╡ 2f62d450-87f7-11eb-144b-b7e2644c1d0e
begin
Y = Int.((stats_2[1]./2).-1) 
X = stats_2[2]
@model logreg(data, X) = begin
    n = 200
	m = 1
    σ ~ truncated(Cauchy(0, 1), 0., Inf)
	#intercept ~ Normal(0, σ)
    β ~ MvNormal(zeros(m), σ)
    Z = logistic.(X.*β) #intercept .+ 
    for i=1:n
        data[i] ~ Bernoulli(Z[i])
    end 
end
chain = sample(logreg(Y, X), NUTS(), 500) 	
plot(chain)
β = hcat(get(chain, :β).β...)
βmean = vec(mean(β, dims=1))
end

# ╔═╡ ea4c2620-8812-11eb-347b-79612e707769
begin
logistic.(X.*βmean)
end

# ╔═╡ ff0cb820-7ff9-11eb-00e7-cf90fe869537
begin
p7 = plot(stats_2[2], stats_2[3], grid=false, color=:black, label="Pop size after t generations")
xlabel!("\$u\$")
ylabel!("Number of individuals")
end

# ╔═╡ f0073e40-7c70-11eb-028a-978df2476189
#Allocates some random probabalities for unreduced gamete formation for 2n, 3n, 4n individuals respectively, for example [0.9,0.1,0.,0.] -> 2n ind has probability of 0.9 to produce 1n gametes, 0.1 for 2n (unreduced) gamete,...this doesn't take into account aneuploid gametes yet.
UG_2 = UnreducedGamete([[0.6,0.4,0.,0.],[0.,0.,0.,0.],[0.,.4,0.0,0.6]])

# ╔═╡ 90771400-7c6b-11eb-33da-4f1bf08039e4
#On the rows/columns (symmetric matrix) are the ploidy levels of the gametes to be combined (1n to 4n).The upper bound is now at 4n, all ploidy levels above are not viable (ie. combination of 2n and 3n gamete has 0 viability).
OV_2 = OffspringViability([[1.,0.,.1,0.],[0.,1.,0.,0.],[.1,0.,0.,0.],[0.,0.,0.,0.]])

# ╔═╡ f89b70d0-7c70-11eb-2d9c-d9abf241dd17
sim_ploidyvar2 = evolving_deme_ploidyvar(d_p,500,UG_2,OV_2)

# ╔═╡ 28e30f50-7c71-11eb-072c-dffe12b4b437
begin
	pf2_p2 = sim_ploidyvar2.p2
	pf3_p2 = sim_ploidyvar2.p3
	pf4_p2 = sim_ploidyvar2.p4
	plot(pf2_p2, grid=false, color=:blue, label="diploid", title="u=v=0.1")
	#plot!(pf3_p2, grid=false, color=:green, label="triploid")
	plot!(pf4_p2, grid=false, color=:red, label="tetraploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ ed0024e0-7cd5-11eb-1548-f70db6fdda06
function grid_search_2(t)
	OV = OffspringViability([[1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,0.]])
	ploidy = []
	param = []
	pop_size = []
	for u in range(0, stop=0.5, length=t)
		UG = UnreducedGamete([[1-u,u,0.,0.],[0.,0.,0.,0.],[0.0,1-u,0.0,u]])
		sim_ploidyvar = evolving_deme_ploidyvar(d_p,50,UG,OV)
		if sim_ploidyvar.p2[end] >= sim_ploidyvar.p4[end]
			push!(ploidy,2)
		else
			push!(ploidy,4)
		end
		pop = (sim_ploidyvar.p2[end] + sim_ploidyvar.p4[end])
		push!(pop_size,pop)
		push!(param, u)
	end
	ploidy, param, pop_size
end

# ╔═╡ 52fc8ae0-7cd6-11eb-1731-f3c71b08b7a9
stats_3 = grid_search_2(500)

# ╔═╡ c88e13d0-7c7e-11eb-0119-731f7ddc65ce
dp_2 = [(stats_2[2][i],stats_2[1][i],stats_3[3][i]) for i in 1:200]

# ╔═╡ 37231930-7c7f-11eb-2524-9115fc85d926
begin
p3 = scatter(dp_2, label=false, title="UG for 2n only")
vline!([0.17], label="u=0.17",linewidth=5)
xlabel!("u (unreduced gamete formation 2n)")
ylabel!("Ploidy")
end

# ╔═╡ 41392210-7c85-11eb-19bf-675f7a929063
binz = bin(dp_2)

# ╔═╡ b7e83d6e-7c84-11eb-0e72-7b4c25184221
begin
	p1 = histogram(binz[1], bins=50, label=false, title="Diploid > Tetraploid")
	vline!([0.17], label="u=0.17",linewidth=5)
	xlabel!("u (unreduced gamete formation 2n)")
	ylabel!("Count")
end

# ╔═╡ 616389e0-7c85-11eb-13c9-fbf0587d3ccc
begin
	p2 = histogram(binz[2], bins=50, label=false, title="Tetraploid > Diploid")
	vline!([0.17], label="u=0.17",linewidth=5)
	xlabel!("u (unreduced gamete formation 2n)")
	ylabel!("Count")
end

# ╔═╡ 7f182630-7cc6-11eb-071d-0d22322e64e2
plot(p1,p2,p3,p7)

# ╔═╡ 99dcebd2-7cd6-11eb-3b96-85bc3a0940d8
dp_3 = [(stats_3[2][i],stats_3[1][i],stats_3[3][i]) for i in 1:500]

# ╔═╡ c0370950-7cd6-11eb-1113-9b2c9104da16
begin
p4 = scatter(dp_3, label=false, title="UG for both 2n and 4n")
vline!([0.2], label="u=0.2",linewidth=5)
xlabel!("u=v (unreduced gamete formation 2n and 4n)")
ylabel!("Ploidy")
end

# ╔═╡ fb49b790-7cd6-11eb-3043-3563217174e9
binz2 = bin(dp_3)

# ╔═╡ 08b7378e-7cd7-11eb-0120-edc0d8f9dbb8
begin
	p5 = histogram(binz2[1], bins=50, label=false, title="Diploid > Tetraploid")
	vline!([0.2], label="u=0.2",linewidth=5)
	xlabel!("u=v")
	ylabel!("Count")
end

# ╔═╡ 38ef6cc0-7cd7-11eb-26e6-1767c49d50de
begin
	p6 = histogram(binz2[2], bins=50, label=false, title="Tetraploid > Diploid")
	vline!([0.2], label="u=0.2",linewidth=5)
	xlabel!("u=v")
	ylabel!("Count")
end

# ╔═╡ 55373330-7ff8-11eb-351c-9df644942e4d
begin
p8 = plot(stats_3[2], stats_3[3], grid=false, color=:black, label="Pop size after t generations")
xlabel!("\$u\$")
ylabel!("Number of individuals")
end

# ╔═╡ c5025790-7cd7-11eb-24f8-51d85dc4b968
plot(p5,p6,p4,p8)

# ╔═╡ Cell order:
# ╟─bd2ab750-7c60-11eb-3cc5-6d57dd448dfd
# ╠═98180cf0-7c6b-11eb-0da9-c55b3cf2a066
# ╟─2dadbbf0-7c5f-11eb-3d09-1fcf2b82ead9
# ╠═472102d2-84a8-11eb-2d6d-eb101b7a707b
# ╟─7df71be0-7c70-11eb-3b92-91457ecf38fd
# ╟─c8734230-7c65-11eb-30bd-4314f44ee3ba
# ╠═8f919600-7c6b-11eb-33b7-bd5a762fb3eb
# ╟─e0789de0-7c6e-11eb-0900-3387c904ea0b
# ╟─e5916552-7c6e-11eb-1868-a7b856e0ad1f
# ╠═1d6f4e60-7c6f-11eb-0be9-9341661ea1d2
# ╠═20380b50-7c6f-11eb-14c9-fb3fc36e318b
# ╠═8caaa5de-7c6f-11eb-387a-05cc2add5483
# ╠═13d994e0-7c70-11eb-1617-c37054ccc63b
# ╟─4cde5ba0-7c79-11eb-23e1-d1ac938c0e6b
# ╟─1f6a4e72-7d48-11eb-0933-eb25e7a5b568
# ╟─e315b920-7c7d-11eb-0f34-5782accf7286
# ╟─39c92860-7c7e-11eb-1efe-13bf99e99d1e
# ╠═5c9c61c0-87f2-11eb-107c-eb6817ebc523
# ╠═ec061ade-8801-11eb-2449-f3788b4ffa9b
# ╠═2f62d450-87f7-11eb-144b-b7e2644c1d0e
# ╠═ea4c2620-8812-11eb-347b-79612e707769
# ╠═c88e13d0-7c7e-11eb-0119-731f7ddc65ce
# ╠═37231930-7c7f-11eb-2524-9115fc85d926
# ╠═d75aee00-7c84-11eb-2d5a-b3ca3c4bd625
# ╟─41392210-7c85-11eb-19bf-675f7a929063
# ╟─b7e83d6e-7c84-11eb-0e72-7b4c25184221
# ╟─616389e0-7c85-11eb-13c9-fbf0587d3ccc
# ╟─ff0cb820-7ff9-11eb-00e7-cf90fe869537
# ╟─7f182630-7cc6-11eb-071d-0d22322e64e2
# ╟─a478224e-7c70-11eb-133b-03e027bf81c8
# ╟─57d84bf0-7c6b-11eb-0e68-f99c717f46e4
# ╠═f0073e40-7c70-11eb-028a-978df2476189
# ╠═90771400-7c6b-11eb-33da-4f1bf08039e4
# ╠═f89b70d0-7c70-11eb-2d9c-d9abf241dd17
# ╟─28e30f50-7c71-11eb-072c-dffe12b4b437
# ╠═ed0024e0-7cd5-11eb-1548-f70db6fdda06
# ╠═52fc8ae0-7cd6-11eb-1731-f3c71b08b7a9
# ╟─99dcebd2-7cd6-11eb-3b96-85bc3a0940d8
# ╠═c0370950-7cd6-11eb-1113-9b2c9104da16
# ╟─55373330-7ff8-11eb-351c-9df644942e4d
# ╠═fb49b790-7cd6-11eb-3043-3563217174e9
# ╟─08b7378e-7cd7-11eb-0120-edc0d8f9dbb8
# ╟─38ef6cc0-7cd7-11eb-26e6-1767c49d50de
# ╟─c5025790-7cd7-11eb-24f8-51d85dc4b968
# ╟─1abbabe0-87fb-11eb-271e-53e02f1ec547
# ╟─86be28f0-84ad-11eb-2ccf-6f60ad079d75
# ╟─7f32ce70-85ab-11eb-1c30-1ffd857a2545
# ╟─f641bd10-85aa-11eb-03fc-0de54af93d35
# ╟─8a935580-85ad-11eb-1c01-fbfe6528c0e2
# ╟─495309e0-85ac-11eb-3231-b3dce7ae61c5
# ╟─aab09870-84b1-11eb-11c6-ffa1224114d8
# ╟─302363ee-7c6c-11eb-0d7b-effeadcdeb60
# ╠═cc1dcf70-7c6c-11eb-3664-b7e70e49723a
# ╠═e419ec50-87fa-11eb-230a-8f6396f0e459
# ╠═0144f6d0-87fb-11eb-0fa7-7149bfb316db
# ╟─41a5cff0-7c6c-11eb-2b03-d7337d430df4
