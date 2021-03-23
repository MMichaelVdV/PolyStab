### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ c0283090-882c-11eb-0f37-c546cb7cb8b9
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes

# ╔═╡ cf15a560-8822-11eb-0789-ad5bb86f481b
md""" Some testing to transfer from array to matrix"""

# ╔═╡ 0d230030-8a4b-11eb-2013-5bb2f8014d53
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

# ╔═╡ 69ffbcf0-89d1-11eb-0e06-4da00973549b


# ╔═╡ a04b3210-89c4-11eb-3c4f-4d2665a7442f
md""" ### Functions"""

# ╔═╡ 84c7bd10-89c4-11eb-2448-4ff2943f6a01
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
d_p = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 0.1, 250, [2], 175, d=2.),randagent_p(0.5, 0.1, 250, [4], 0, d=4.)), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0.9 0.1 0. 0. ; 0. 0. 0. 0. ; 0. 0.9 0. 0.1] )

# ╔═╡ ebcf3950-8836-11eb-1092-9fe734282054
d_p.agents

# ╔═╡ eb74ab50-8838-11eb-25e5-010384ccca19
d_p.UG[3,:]

# ╔═╡ ab91dec0-8835-11eb-2594-bb2fed4f14e7
a = randagent_p(0.5, 0.25, 100, 4, d=4.)

# ╔═╡ b86086ae-8835-11eb-225b-bf9c29114e21
size(a.loci)

# ╔═╡ 58da38c0-8836-11eb-18c7-e3087770b4fb
length(a)

# ╔═╡ c127ec50-8837-11eb-21d1-0b4c07df9464
sum(a)

# ╔═╡ d21432d0-8837-11eb-11bc-49a4315c59ff
sum(a.loci)

# ╔═╡ cab6126e-883b-11eb-2fc4-65f85a389ffd
a.loci[1,:]

# ╔═╡ 3c6ace60-883c-11eb-215f-f95e755463d9
a.loci[shuffle(1:end), :]

# ╔═╡ 72920d0e-8908-11eb-1fc7-7f2ec0184880
a1 = randagent_p(0.5, 0.25, 100, 1, d=1.)

# ╔═╡ 72a93e90-8908-11eb-37f8-cd6ad67826a4
a2 = randagent_p(0.5, 0.25, 100, 2, d=2.)

# ╔═╡ c636e860-883a-11eb-27e7-c9e3f70fbe75
num = sample([1.,2.,3.,4.],weights(d_p.UG[ploidy(a)-1,:]))

# ╔═╡ 324f5502-883b-11eb-15de-4fd9dbdbe959
loci = (Float64[], num)

# ╔═╡ 8542bb30-8a4d-11eb-0189-0df3acaa1bae
begin
	
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
	
end

# ╔═╡ 34fecc50-8835-11eb-0e1f-4b2e13ff740e
begin
	
expected_heterozygosity(H₀, t, N) = ((1.0-1.0/N)^t)*H₀
	
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
        		if ag.loci[chr, loc] != 0
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
	map(p->2*p*(1-p), freqs)
end

"""
	evolving_deme_ploidyvar(d::AbstractDeme, ngen)
Simulate a single deme with mixed ploidy, malthusian fitness and unreduced gamete formation.
"""
function evolving_deme_ploidyvar(d::MixedPloidyDeme, ngen; heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = allelefreqs_p, pf = ploidy_freq)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[1]]
	p3 = [ploidy_freq(d)[2]]
	p4 = [ploidy_freq(d)[3]]
	
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
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, het=het,tm=tm, af=af) 
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

end

# ╔═╡ 4fd07fa0-8908-11eb-11aa-c5f5c06cec55
viability(a1, a2, d_p) 

# ╔═╡ c6c83340-8837-11eb-048c-f91ba08f9826
trait(a)

# ╔═╡ 0e650dc0-8835-11eb-122e-f1b02f2471a0
sim_ploidyvar1 = evolving_deme_ploidyvar(d_p,500)

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

# ╔═╡ 727de160-89cf-11eb-2111-bf2ea7fef8b7
begin
	traitmean_ploidy = map(mean, sim_ploidyvar1.tm)
	plot(traitmean_ploidy, grid=false, color=:black, label=false, title="Effect of stabilizing selection")
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ a260fbd0-8a4a-11eb-03a6-0d789a020d9b
function grid_search(t)
	ploidy = []
	param = []
	pop_size = []
	for u in range(0, stop=0.5, length=t)
		d_p = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 0.25, 100, [2], 50, d=2.),randagent_p(0.5, 0.25, 100, [4], 0, d=4.)), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.] )
		sim_ploidyvar = evolving_deme_UG(d_p,50)
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

# ╔═╡ bd382b40-8a4a-11eb-2069-51fd29792863
stats_2 = grid_search(200)

# ╔═╡ dc085a40-8a4a-11eb-3c18-8fb23db823fd
dp_2 = [(stats_2[2][i],stats_2[1][i],stats_2[3][i]) for i in 1:200]

# ╔═╡ de14f320-8a4a-11eb-3fa1-df9b2d93a4b5
begin
p3 = scatter(dp_2, label=false, title="UG for 2n only")
vline!([0.17], label="u=0.17",linewidth=5)
xlabel!("u (unreduced gamete formation 2n)")
ylabel!("Ploidy")
end

# ╔═╡ 3a846eb0-8a4b-11eb-2c21-b772a6f5a9ae
binz = bin(dp_2)

# ╔═╡ 0f8421b0-8a4b-11eb-32d8-ab006d93bcdb
begin
	p12 = histogram(binz[1], bins=50, label=false, title="Diploid > Tetraploid")
	vline!([0.17], label="u=0.17",linewidth=5)
	xlabel!("u (unreduced gamete formation 2n)")
	ylabel!("Count")
end

# ╔═╡ 1ab14ea0-8a4b-11eb-3d46-db34e82fd19c
begin
	p13 = histogram(binz[2], bins=50, label=false, title="Tetraploid > Diploid")
	vline!([0.17], label="u=0.17",linewidth=5)
	xlabel!("u (unreduced gamete formation 2n)")
	ylabel!("Count")
end

# ╔═╡ 41635ab0-8a4c-11eb-3857-e7465c7028d7
begin
p14 = scatter(dp_2, label=false, title="UG for 2n only")
vline!([0.17], label="u=0.17",linewidth=5)
xlabel!("u (unreduced gamete formation 2n)")
ylabel!("Ploidy")
end

# ╔═╡ c3953840-8836-11eb-2b8b-516d9074e45d
ploidy_freq(d_p)

# ╔═╡ d15eb34e-88d8-11eb-145a-0f40e5907c42
recombine_poly(d_p)

# ╔═╡ f6e3deb0-8906-11eb-273d-b127626fa404
mate_p(a,a,d_p)

# ╔═╡ b9bf2800-89c4-11eb-37a9-bbdf6bf60848
begin
	Hₒ_ploidyvar = map(mean, sim_ploidyvar1.het)
	plot(Hₒ_ploidyvar, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH mixed ploidy with selection")
	plot!(0:sim_ploidyvar1.ngen+1, 
		t->expected_heterozygosity(2*0.5*(0.5), t, 100),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 0fae82f0-89d0-11eb-302e-e557072f42df
allelefreqs_p(d_p)

# ╔═╡ 135802f0-89d0-11eb-2332-49218f26b70d
heterozygosities_p(d_p)

# ╔═╡ b6652790-89c4-11eb-0e6c-afb101216229
begin
	
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

end

# ╔═╡ 04fd6e62-89d1-11eb-0574-216f9a2482cb
begin
	
habi = Habitat(demes=[d_p])	
g_lin = linear_gradient(habi)	
hab = initiate_habitat(g_lin, d_p, 0.5, 0.1, 250, 50)	
sim_hab = evolving_habitat(hab, 100)	

	
end


# ╔═╡ 3e9da0e0-89d1-11eb-1614-699289c788b3
begin 

K = d_p.K
σ = hab.σ
b = hab.b
Vs = d_p.Vs
rm = d_p.rm
Dm = hab.Dm
s = 1
pop_sizes = [length(deme) for deme  in sim_hab[1].demes]
ppf1 = [ploidy_freq(deme)[1] for deme  in sim_hab[1].demes]
ppf2 = [ploidy_freq(deme)[2] for deme  in sim_hab[1].demes]
ppf3 = [ploidy_freq(deme)[3] for deme  in sim_hab[1].demes]
margin = (sqrt(2)*b*σ)/((2*rm*sqrt(Vs))-b*σ) .>= 0.15.*pop_sizes.*σ*sqrt(s)
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
trait_means = [trait_mean(deme) for deme in sim_hab[1].demes]
	
function f_het_demes(h::Habitat)
	het_demes = []
	cordsh = []
	for (i, deme) in enumerate(h.demes)
	if length(deme) != 0
		het = 0.5^2*mean(heterozygosities_p(deme)) #why 0.5^2 ??
		push!(cordsh,i)
		push!(het_demes,het)
	else
		push!(cordsh,i)
		push!(het_demes,0)
		end
	end
	het_demes, cordsh
end

end

# ╔═╡ 9d547460-89d1-11eb-0993-e33e38342212
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

# ╔═╡ a4b18680-89d1-11eb-38e5-333f4851ea3f
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
	gif(anim_range, "fizzypopM.gif", fps = 3)
end

# ╔═╡ e496a460-89d1-11eb-11a4-a9d35cf68963
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
		#plot!(ppf1, grid=false, color=:blue, label="Diploid")
		#plot!(ppf2, grid=false, color=:green, label="Triploid")
		#plot!(ppf3, grid=false, color=:red, label="Tetraploid")
		hline!([b*σ*sqrt(Vs)], label = "E(V_G)")
		
		
		xlabel!("Space")
		ylabel!("Vg")
	end every 1
	gif(anim_range_Vg, "fizzyVgM.gif", fps = 3)
end

# ╔═╡ ee7928e0-89d1-11eb-38b2-df6427fb0853
begin
	anim_range_trait = @animate for i ∈ 1:100
		sim_habA = sim_hab.data[i]
		trait_means = [trait_mean(deme) for deme in sim_habA.demes]
	    trait_means_p = map(mean, trait_means)
		trait_agents, cordst = f_trait_agents(sim_habA)
		#p2 = plot(trait_means_p, grid=false, color=:black, label="Z_mean deme")
		plot(g_lin, grid=false, color=:blue, label="Z optimum", linestyle=:dash)
		plot!(cordst,trait_agents, label="Z agents")
		xlabel!("Space")
		ylabel!("Trait Z")
	end every 1
	gif(anim_range_trait, "fizzyM.gif", fps = 3)
end

# ╔═╡ Cell order:
# ╟─cf15a560-8822-11eb-0789-ad5bb86f481b
# ╠═c0283090-882c-11eb-0f37-c546cb7cb8b9
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
# ╠═a260fbd0-8a4a-11eb-03a6-0d789a020d9b
# ╠═bd382b40-8a4a-11eb-2069-51fd29792863
# ╠═dc085a40-8a4a-11eb-3c18-8fb23db823fd
# ╠═de14f320-8a4a-11eb-3fa1-df9b2d93a4b5
# ╠═0d230030-8a4b-11eb-2013-5bb2f8014d53
# ╠═3a846eb0-8a4b-11eb-2c21-b772a6f5a9ae
# ╠═0f8421b0-8a4b-11eb-32d8-ab006d93bcdb
# ╠═1ab14ea0-8a4b-11eb-3d46-db34e82fd19c
# ╠═41635ab0-8a4c-11eb-3857-e7465c7028d7
# ╠═c3953840-8836-11eb-2b8b-516d9074e45d
# ╠═eb74ab50-8838-11eb-25e5-010384ccca19
# ╠═c636e860-883a-11eb-27e7-c9e3f70fbe75
# ╠═324f5502-883b-11eb-15de-4fd9dbdbe959
# ╠═cab6126e-883b-11eb-2fc4-65f85a389ffd
# ╠═3c6ace60-883c-11eb-215f-f95e755463d9
# ╠═d15eb34e-88d8-11eb-145a-0f40e5907c42
# ╠═f6e3deb0-8906-11eb-273d-b127626fa404
# ╠═727de160-89cf-11eb-2111-bf2ea7fef8b7
# ╠═b9bf2800-89c4-11eb-37a9-bbdf6bf60848
# ╠═0fae82f0-89d0-11eb-302e-e557072f42df
# ╠═135802f0-89d0-11eb-2332-49218f26b70d
# ╠═04fd6e62-89d1-11eb-0574-216f9a2482cb
# ╠═69ffbcf0-89d1-11eb-0e06-4da00973549b
# ╠═3e9da0e0-89d1-11eb-1614-699289c788b3
# ╠═9d547460-89d1-11eb-0993-e33e38342212
# ╠═a4b18680-89d1-11eb-38e5-333f4851ea3f
# ╠═e496a460-89d1-11eb-11a4-a9d35cf68963
# ╠═ee7928e0-89d1-11eb-38b2-df6427fb0853
# ╟─a04b3210-89c4-11eb-3c4f-4d2665a7442f
# ╠═84c7bd10-89c4-11eb-2448-4ff2943f6a01
# ╠═34fecc50-8835-11eb-0e1f-4b2e13ff740e
# ╠═b6652790-89c4-11eb-0e6c-afb101216229
# ╠═8542bb30-8a4d-11eb-0189-0df3acaa1bae
