### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ d2a4995c-c62a-47d7-bcc2-626f8989bf3d
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, StatsPlots

# ╔═╡ fb30b92f-2cb8-406a-b50d-238282889028
using PolyStab

# ╔═╡ b7182aa0-1635-4078-b736-083fe0e446fb
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, evolving_neutraldeme, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy, trait_mean, randagent, evolving_haploiddeme, mate_p, evolving_selectiondemeh, malthusian_fitness, number_of_offspring, ismock, ploidy_freq, f_trait_agents, mutate

# ╔═╡ 381a9880-e3cd-11eb-0b70-7bb37ea370b2
using PolyStab: IslandDeme, directional_selection, popsize

# ╔═╡ c5a678de-2382-4ab2-98ef-a91cf66e2d9f
Poisson(0.01)

# ╔═╡ 19e0c7e8-de79-46fd-85f6-822d586be8b9
"""
	mating_PnB_s(d::MixedPloidyDeme{A})

Mating in a mixed ploidy deme with unreduced gamete formation and partner
choice weighted by malthusian fitness.
"""
function mating_PnB_s(d::IslandDeme{A}, s::Float64) where A
	new_agents = A[]
	fitnesses = exp.(directional_selection(d))
	pl = ploidy.(d.agents)
	for i=1:length(d)
		B1 = d.agents[i]
		if rand() < s 
			pls = (pl .*0)
			pls[i] = 1.
		else pls = (pl .*0) .+ 1.
		end
		fits = pls.*fitnesses
		#fitnesses = exp.(malthusian_fitness(d) .+ (-s .* (trait(B1) .- trait.(d.agents)).^2))
		noff = number_of_offspring(d, B1)
		Bs = sample(d.agents, weights(fits), noff) 
        offspring = filter(!ismock, map(B2->mate_p(B1, B2, d), Bs))
        new_agents = vcat(new_agents, offspring)
    end
    d(new_agents)
end

# ╔═╡ 3f31c2ce-231c-492e-9291-e044f7a9ee6c
"""
	mating_PnB_a(d::MixedPloidyDeme{A})

Mating in a mixed ploidy deme with unreduced gamete formation and partner
choice weighted by malthusian fitness.
"""
function mating_PnB_a(d::IslandDeme{A}, a::Float64) where A
	new_agents = A[]
	fitnesses = exp.(directional_selection(d))
	pl = ploidy.(d.agents)
	for i=1:length(d)
		B1 = d.agents[i]
		if rand() < a 
			pls = (pl .*0) .+ (pl.==ploidy(B1))
		else pls = (pl .*0) .+ 1.
		end
		fits = pls.*fitnesses
		#fitnesses = exp.(malthusian_fitness(d) .+ (-s .* (trait(B1) .- trait.(d.agents)).^2))
		noff = popsize(d) < 350 ? number_of_offspring(d,B1) : 1 #limits popsize to prevent exponential explosion
		Bs = sample(d.agents, weights(fits), noff) 
        offspring = filter(!ismock, map(B2->mate_p(B1, B2, d), Bs))
        new_agents = vcat(new_agents, offspring)
    end
    d(new_agents)
end

# ╔═╡ aa650c14-ff9f-4f6b-935e-7c9ecf2ff823
function evolving_islandwbreak_a(d::IslandDeme, M, L, s, ngen; 
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
	c = 0
	
	for t=1:ngen
		if popsize(d) < L
		migrants = rand(Poisson(M))
		for m in 1:migrants
			migrant = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1)[1]
			push!(d.agents,migrant)
		end
		d = mating_PnB_a(d,s)
		#d = mutate(d) 
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[2])
		push!(p3, ploidy_freq(d)[3])
		push!(p4, ploidy_freq(d)[4])
		push!(fta, f_trait_agents(d))
		c += 1
		else break
		end
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=c, het=het,tm=tm, af=af, fta=fta) 
end

# ╔═╡ 9d92f46e-a780-4e64-9846-4ab761cd84d3
us = 0.2

# ╔═╡ 8e760999-921a-41ef-8114-97501c0e46e7
UGs = [0. 0. 0. 0. ; 1-us us 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.]

# ╔═╡ d552fb3d-0107-4a6e-b6e4-bdf966e81d87
a = .5

# ╔═╡ 475064f2-b2f1-4dba-976b-d4b75cf40547
#dz1, m=[0.01,0.1,1,10], u=0.01
begin
est_1_001_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=13.5, β=0.25), 0.01, 100, a, 10000)), 1:20)
est_1_01_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=13.5, β=0.25), 0.1, 100, a, 10000)), 1:20)
est_1_1_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=13.5, β=0.25), 1., 100, a, 10000)), 1:20)
est_1_10_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=13.5, β=0.25), 10., 100, a, 10000)), 1:20)
end

# ╔═╡ 03af42af-5ce4-4a97-8270-5af2338f351a
#dz2, m=[0.01,0.1,1,10], u=0.01
begin
est_2_001_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=14.5, β=0.25), 0.01, 100, a, 10000)), 1:20)
est_2_01_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=14.5, β=0.25), 0.1, 100, a, 10000)), 1:20)
est_2_1_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG =UGs, μ=0., θ=14.5, β=0.25), 1., 100, a, 10000)), 1:20)
est_2_10_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=14.5, β=0.25), 10., 100, a, 10000)), 1:20)
end

# ╔═╡ f80deba9-6f73-4a39-96b4-8b4fcb587065
#dz3, m=[0.01,0.1,1,10], u=0.01
begin
est_3_001_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=15.5, β=0.25), 0.01, 100, a, 10000)), 1:20)
est_3_01_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=15.5, β=0.25), 0.1, 100, a, 10000)), 1:20)
est_3_1_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=15.5, β=0.25), 1., 100, a, 10000)), 1:20)
est_3_10_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=15.5, β=0.25), 10., 100, a, 10000)), 1:20)
end

# ╔═╡ 8bb17e05-c7a0-405b-861d-0b8f1dc01045
#dz4, m=[0.01,0.1,1,10], u=0.01
begin
est_4_001_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=16.5, β=0.25), 0.01, 100, a, 10000)), 1:20)
est_4_01_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=16.5, β=0.25), 0.1, 100, a, 10000)), 1:20)
est_4_1_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=16.5, β=0.25), 1., 100, a, 10000)), 1:20)
est_4_10_001 = map(x->(evolving_islandwbreak_a(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=16.5, β=0.25), 10., 100, a, 10000)), 1:20)
end

# ╔═╡ 3785c38d-ea09-435d-ad5f-da86e5e6f4dc
begin
	p1_001=sum([est_1_001_001[i].p4[end] for i in 1:20] .> 100)/sum([est_1_001_001[i].pop[end] for i in 1:20] .> 100)
	p1_01=sum([est_1_01_001[i].p4[end] for i in 1:20] .> 100)/sum([est_1_01_001[i].pop[end] for i in 1:20] .> 100)
	p1_1=sum([est_1_1_001[i].p4[end] for i in 1:20] .> 100)/sum([est_1_1_001[i].pop[end] for i in 1:20] .> 100)
	p1_10=sum([est_1_10_001[i].p4[end] for i in 1:20] .> 100)/sum([est_1_10_001[i].pop[end] for i in 1:20] .> 100)
	
	p2_001=sum([est_2_001_001[i].p4[end] for i in 1:20] .> 100)/sum([est_2_001_001[i].pop[end] for i in 1:20] .> 100)
	p2_01=sum([est_2_01_001[i].p4[end] for i in 1:20] .> 100)/sum([est_2_01_001[i].pop[end] for i in 1:20] .> 100)
	p2_1=sum([est_2_1_001[i].p4[end] for i in 1:20] .> 100)/sum([est_2_1_001[i].pop[end] for i in 1:20] .> 100)
	p2_10=sum([est_2_10_001[i].p4[end] for i in 1:20] .> 100)/sum([est_2_10_001[i].pop[end] for i in 1:20] .> 100)
	
	p3_001=sum([est_3_001_001[i].p4[end] for i in 1:20] .> 100)/sum([est_3_001_001[i].pop[end] for i in 1:20] .> 100)
	p3_01=sum([est_3_01_001[i].p4[end] for i in 1:20] .> 100)/sum([est_3_01_001[i].pop[end] for i in 1:20] .> 100)
	p3_1=sum([est_3_1_001[i].p4[end] for i in 1:20] .> 100)/sum([est_3_1_001[i].pop[end] for i in 1:20] .> 100)
	p3_10=sum([est_3_10_001[i].p4[end] for i in 1:20] .> 100)/sum([est_3_10_001[i].pop[end] for i in 1:20] .> 100)
	
	p4_001=sum([est_4_001_001[i].p4[end] for i in 1:20] .> 100)/sum([est_4_001_001[i].pop[end] for i in 1:20] .> 100)
	p4_01=sum([est_4_01_001[i].p4[end] for i in 1:20] .> 100)/sum([est_4_01_001[i].pop[end] for i in 1:20] .> 100)
	p4_1=sum([est_4_1_001[i].p4[end] for i in 1:20] .> 100)/sum([est_4_1_001[i].pop[end] for i in 1:20] .> 100)
	p4_10=sum([est_4_10_001[i].p4[end] for i in 1:20] .> 100)/sum([est_4_10_001[i].pop[end] for i in 1:20] .> 100)
end

# ╔═╡ a0c270c9-6a85-4137-9bc5-4c32bcbbf876
begin
	sep1_001 = sqrt(p1_001*(1-p1_001))/20
	sep1_01 = sqrt(p1_01*(1-p1_01))/20
	sep1_1 = sqrt(p1_1*(1-p1_1))/20
	sep1_10 = sqrt(p1_10*(1-p1_10))/20
	
	sep2_001 = sqrt(p2_001*(1-p2_001))/20
	sep2_01 = sqrt(p2_01*(1-p2_01))/20
	sep2_1 = sqrt(p2_1*(1-p2_1))/20
	sep2_10 = sqrt(p2_10*(1-p2_10))/20
	
	sep3_001 = sqrt(p3_001*(1-p3_001))/20
	sep3_01 = sqrt(p3_01*(1-p3_01))/20
	sep3_1 = sqrt(p3_1*(1-p3_1))/20
	sep3_10 = sqrt(p3_10*(1-p3_10))/20
	
	sep4_001 = sqrt(p4_001*(1-p4_001))/20
	sep4_01 = sqrt(p4_01*(1-p4_01))/20
	sep4_1 = sqrt(p4_1*(1-p4_1))/20
	sep4_10 = sqrt(p4_10*(1-p4_10))/20
end

# ╔═╡ 02c3be6a-32f8-4cca-8397-b779899938bb
begin
estab2_01 = plot([1,2,3,4],[p1_001,p1_01,p1_1,p1_10], label="△z=1", marker = ([:hex :d]), color=:green, title="u=0.2, a=0.5", yerror=[sep1_001, sep1_01, sep1_1, sep1_10])
plot!([1,2,3,4],[p2_001,p2_01,p2_1,p2_10], label="△z=2", marker = ([:hex :d]), color=:blue, yerror=[sep2_001, sep2_01, sep2_1, sep2_10])
plot!([1,2,3,4],[p3_001,p3_01,p3_1,p3_10], label="△z=3", marker = ([:hex :d]), color=:red, yerror=[sep3_001, sep3_01, sep3_1, sep3_10])
plot!([1,2,3,4],[p4_001,p4_01,p4_1,p4_10], label="△z=4", marker = ([:hex :d]), color=:black, yerror=[sep4_001, sep4_01, sep4_1, sep4_10])
#plot!([1,2,3,4],[p5_001,p5_01,p5_1,p5_10], label="△z=5", marker = ([:hex :d]), color=:blue, ribbon=[sep5_001, sep5_01, sep5_1, sep5_10])
#plot!([1,2,3,4],[p6_001,p6_01,p6_1,p6_10], label="△z=6", marker = ([:hex :d]), color=:purple, ribbon=[sep6_001, sep6_01, sep6_1, sep6_10])
xlabel!("Migration rate")
ylabel!("Probability of tetraploid establishment")
#savefig(estab_01, "Estab_01")
end

# ╔═╡ 99abbd3d-5b5b-4482-bd41-afe433cf1ae7
begin
st1_001=mean([est_1_001_001[i].ngen for i in 1:20])
st1_01=mean([est_1_01_001[i].ngen for i in 1:20])
st1_1=mean([est_1_1_001[i].ngen for i in 1:20])
st1_10=mean([est_1_10_001[i].ngen for i in 1:20])
se1_001=sem([est_1_001_001[i].ngen for i in 1:20])
se1_01=sem([est_1_01_001[i].ngen for i in 1:20])
se1_1=sem([est_1_1_001[i].ngen for i in 1:20])
se1_10=sem([est_1_10_001[i].ngen for i in 1:20])

st2_001=mean([est_2_001_001[i].ngen for i in 1:20])
st2_01=mean([est_2_01_001[i].ngen for i in 1:20])
st2_1=mean([est_2_1_001[i].ngen for i in 1:20])
st2_10=mean([est_2_10_001[i].ngen for i in 1:20])
se2_001=sem([est_2_001_001[i].ngen for i in 1:20])
se2_01=sem([est_2_01_001[i].ngen for i in 1:20])
se2_1=sem([est_2_1_001[i].ngen for i in 1:20])
se2_10=sem([est_2_10_001[i].ngen for i in 1:20])

st3_001=mean([est_3_001_001[i].ngen for i in 1:20])
st3_01=mean([est_3_01_001[i].ngen for i in 1:20])
st3_1=mean([est_3_1_001[i].ngen for i in 1:20])
st3_10=mean([est_3_10_001[i].ngen for i in 1:20])
se3_001=sem([est_3_001_001[i].ngen for i in 1:20])
se3_01=sem([est_3_01_001[i].ngen for i in 1:20])
se3_1=sem([est_3_1_001[i].ngen for i in 1:20])
se3_10=sem([est_3_10_001[i].ngen for i in 1:20])

st4_001=mean([est_4_001_001[i].ngen for i in 1:20])
st4_01=mean([est_4_01_001[i].ngen for i in 1:20])
st4_1=mean([est_4_1_001[i].ngen for i in 1:20])
st4_10=mean([est_4_10_001[i].ngen for i in 1:20])
se4_001=sem([est_4_001_001[i].ngen for i in 1:20])
se4_01=sem([est_4_01_001[i].ngen for i in 1:20])
se4_1=sem([est_4_1_001[i].ngen for i in 1:20])
se4_10=sem([est_4_10_001[i].ngen for i in 1:20])
end

# ╔═╡ abbca9c1-f65b-4f88-9ebb-713436e70314
begin
estab1_01 = plot([1,2,3,4],[st1_001, st1_01, st1_1, st1_10], label="△z=1", marker = ([:hex :d]), color=:green, title="u=0.2, a=.5", yerror=[se1_001, se1_01, se1_1, se1_10])
plot!([1,2,3,4],[st2_001, st2_01, st2_1, st2_10], label="△z=2", marker = ([:hex :d]), color=:blue, yerror=[se2_001, se2_01, se2_1, se2_10])
plot!([1,2,3,4],[st3_001, st3_01, st3_1, st3_10], label="△z=3", marker = ([:hex :d]), color=:red, yerror=[se3_001, se3_01, se3_1, se3_10])
plot!([1,2,3,4],[st4_001, st4_01, st4_1, st4_10], label="△z=4", marker = ([:hex :d]), color=:black, yerror=[se4_001, se4_01, se4_1, se4_10])
#plot!([1,2,3,4],[st5_001, st5_01, st5_1, st5_10], label="△z=5", marker = ([:hex :d]), color=:orange, ribbon=[se5_001, se5_01, se5_1, se5_10])
#plot!([1,2,3,4],[st6_001, st6_01, st6_1, st6_10], label="△z=6", marker = ([:hex :d]), color=:red, ribbon=[se6_001, se6_01, se6_1, se6_10])
xlabel!("Migration rate")
ylabel!("Time to establish")
#savefig(p2nsimUG, "Estab_2n_UG020020")
end

# ╔═╡ 331c6874-a968-46d1-87aa-32a6019ceb5b
begin
	c = plot(estab1_01, estab2_01)
	#savefig(c, "Estab_005")
end

# ╔═╡ Cell order:
# ╠═d2a4995c-c62a-47d7-bcc2-626f8989bf3d
# ╠═c5a678de-2382-4ab2-98ef-a91cf66e2d9f
# ╠═fb30b92f-2cb8-406a-b50d-238282889028
# ╠═b7182aa0-1635-4078-b736-083fe0e446fb
# ╠═381a9880-e3cd-11eb-0b70-7bb37ea370b2
# ╠═19e0c7e8-de79-46fd-85f6-822d586be8b9
# ╠═3f31c2ce-231c-492e-9291-e044f7a9ee6c
# ╠═aa650c14-ff9f-4f6b-935e-7c9ecf2ff823
# ╠═8e760999-921a-41ef-8114-97501c0e46e7
# ╠═9d92f46e-a780-4e64-9846-4ab761cd84d3
# ╠═d552fb3d-0107-4a6e-b6e4-bdf966e81d87
# ╠═475064f2-b2f1-4dba-976b-d4b75cf40547
# ╠═03af42af-5ce4-4a97-8270-5af2338f351a
# ╠═f80deba9-6f73-4a39-96b4-8b4fcb587065
# ╠═8bb17e05-c7a0-405b-861d-0b8f1dc01045
# ╠═3785c38d-ea09-435d-ad5f-da86e5e6f4dc
# ╠═a0c270c9-6a85-4137-9bc5-4c32bcbbf876
# ╠═02c3be6a-32f8-4cca-8397-b779899938bb
# ╠═99abbd3d-5b5b-4482-bd41-afe433cf1ae7
# ╠═abbca9c1-f65b-4f88-9ebb-713436e70314
# ╠═331c6874-a968-46d1-87aa-32a6019ceb5b
