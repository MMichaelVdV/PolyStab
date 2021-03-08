"""
    Agent{T}
Should be an agent of arbitrary (allowed) ploidy level.
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
Base.sum(a::Agent) = sum(sum(a.loci))

ploidy(d::AbstractDeme) = [length(a.loci) for a in d.agents]

"""
"""
function mate(a::Agent, b::Agent)
    newloci = similar(a.loci)
    for i in 1:length(a)
		@inbounds newloci[i] = rand() < 0.5 ? a[i] : b[i]
    end
	Agent(newloci, a.d)
end	

"""
"""
function random_mating(d::AbstractDeme)
    newdeme = similar(d.agents)
    for i=1:length(d)
		@inbounds newdeme[i] = mate(rand(d, 2)...)
	end
	d(newdeme)
end 

"""
"""
function random_mating_p(d::AbstractDeme)
    newdeme = similar(d.agents)
    for i=1:length(d)
		@inbounds newdeme[i] = mate_p(rand(d, 2)...)
	end
	d(agents=newdeme,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end

"""
#This mating function assumes that different cytotypes are incompatible (i.e. when #individuals with different ploidy hybridize, they generate no offspring). This is in #accordance with the model from Levy (1975). Selfing is allowed without cost. The #dynamics might be different when you implement mating alowing for sterile hybrid #offspring that can compete for space (and affects the malthusian fitness).
"""	
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

"""
#This mating function assumes that different cytotypes can be compatible with a decrease in viability (cfr. OffspringViability matrix) (i.e. when #individuals with different ploidy hybridize, they generate have a probability p to generate no viable offspring). Selfing is allowed without cost. This influences the dynamics by including hybrid offspring that can compete for space (and affects the malthusian fitness/density dependence of selection). Still need to incorporate a decent recombination method. 
"""	
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

"""
"""
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

"""
"""
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
	
"""
"""
function allelefreqs(d::AbstractDeme)
    # assuming 0 is absence of allele
    f = j->mapreduce(i->d[i][j] != 0., +, 1:length(d))/length(d)
	frequencies = map(f, 1:length(d[1]))
end

"""
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
"""
function heterozygosities(d::AbstractDeme) 
	freqs=allelefreq(d)
    map(p->p*(1-p), freqs)
end

"""
"""
function heterozygosities_p(d::AbstractDeme)
	freqs=allelefreqs_p(d)
    map(p->p*(1-p), freqs)
end

"""
"""
expected_heterozygosity(H₀, t, N) = ((1.0-1.0/N)^t)*H₀


"""
"""	
function trait(a::Agent)
	sum(a)/a.d
end

"""
"""
function trait_mean(d::AbstractDeme)
	z = Float64[]
	for agent in d.agents
		push!(z,trait(agent))
	end
	sum(z)/length(d)
end

"""
"""
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

"""
"""
function unreduced_gamete(d::AbstractDeme{A}) where A
	new_agents =  A[]
	for agent in d.agents
		push!(new_agents, unreduced_gamete(d,agent))
	end
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end	

"""
"""
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

"""
"""
function malthusian_fitness(d::AbstractDeme,a::Agent)
    N = length(d)
    z = trait(a)
    return d.rm*(1-(N/d.K))-((z-d.θ)^2)/(2*d.Vs)
end 

"""
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
"""
number_of_offspring(d::AbstractDeme,a::Agent) = rand(Poisson(exp(malthusian_fitness(d::AbstractDeme,a::Agent))))

"""
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
		if mate_p(B1,B2) != 0
			for c in 1:noff 
				push!(new_agents, mate_p(B1,B2))
			end
		end
	end
	MixedPloidyDeme(agents=new_agents,K=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end

"""
"""
function viability(a::Agent, b::Agent, OV::OffspringViability)
	return OV.viability[ploidy(a)][ploidy(b)]
end

"""
"""
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

"""
"""
function unreduced_gamete(d::AbstractDeme{A}, UG::UnreducedGamete) where A
	new_agents =  A[]
	for agent in d.agents
		push!(new_agents, unreduced_gamete(agent,UG))
	end
	MixedPloidyDeme(agents=new_agents,k=d.K,θ=d.θ,rm=d.rm,Vs=d.Vs,u=d.u,μ=d.μ)
end		

"""
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
	MixedPloidyDeme(agents=new_agents)
end	

"""
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