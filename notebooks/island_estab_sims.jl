### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ f6f9d8d9-e253-4431-9690-802511aa869f
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, StatsPlots

# ╔═╡ d4c6dd9e-b993-49cc-910b-7d408c5fe500
using PolyStab: Agent, randagent_p, MixedPloidyDeme, IslandDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, random_mating, evolving_neutraldeme, evolving_islanddeme, malthusian_fitness, trait_mean, ploidy_freq, f_trait_agents, mating_PnB_x, directional_selection, mating_PnB, mutate, popsize

# ╔═╡ 9aa821e0-a906-11eb-10e2-89580ee28b2c
md""" ### The effect of migration load on polyploid establishment"""

# ╔═╡ e9d83c94-d3c4-41ca-b4d8-30ee1eda413d
us = 0.05

# ╔═╡ fb6945a6-5c05-474d-b8f8-2749fd304a30
UGs = [0. 0. 0. 0. ; 1-us us 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.]

# ╔═╡ 3609e21d-da26-4cb6-93fe-d3632d19d3d7
#dz4, m=[0.01,0.1,1,10], u=0.01
#begin
#est_5_001_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=17.5, β=0.25), 0.01, 100, 10000)), 1:20)
#est_5_01_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=17.5, β=0.25), 0.1, 100, 10000)), 1:20)
#est_5_1_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=17.5, β=0.25), 1., 100, 10000)), 1:20)
#est_5_10_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=17.5, β=0.25), 10., 100, 10000)), 1:20)
#end

# ╔═╡ 95c2aa35-fb67-4dea-9cf0-d7d8d8b4d31c
#dz4, m=[0.01,0.1,1,10], u=0.01
#begin
#est_6_001_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=17.5, β=0.25), 0.01, 100, 10000)), 1:20)
#est_6_01_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=17.5, β=0.25), 0.1, 100, 10000)), 1:20)
#est_6_1_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=17.5, β=0.25), 1., 100, 10000)), 1:20)
#est_6_10_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=17.5, β=0.25), 10., 100, 10000)), 1:20)
#end

# ╔═╡ b690a6ed-5e38-40b4-ac80-e2f9ed3882c4
md""" #### Functions"""

# ╔═╡ adbeed31-db8c-46ce-9c84-ee48d2173c9c
function evolving_islandwbreak(d::IslandDeme, M, L, ngen; 
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
		d = mating_PnB(d)
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

# ╔═╡ 5cfb3522-997e-42c3-8156-26fe6bf615c0
#dz1, m=[0.01,0.1,1,10], u=0.01
begin
est_1_001_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=13.5, β=0.25), 0.01, 100, 10000)), 1:20)
est_1_01_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=13.5, β=0.25), 0.1, 100, 10000)), 1:20)
est_1_1_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=13.5, β=0.25), 1., 100, 10000)), 1:20)
est_1_10_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=13.5, β=0.25), 10., 100, 10000)), 1:20)
end

# ╔═╡ 5e04a741-c863-4689-95f0-e55370e53d7d
[est_1_001_001[i].ngen for i in 1:20]

# ╔═╡ 41ffb48d-ffdf-4b24-bbf6-4c3194295d71
#dz2, m=[0.01,0.1,1,10], u=0.01
begin
est_2_001_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=14.5, β=0.25), 0.01, 100, 10000)), 1:20)
est_2_01_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=14.5, β=0.25), 0.1, 100, 10000)), 1:20)
est_2_1_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG =UGs, μ=0., θ=14.5, β=0.25), 1., 100, 10000)), 1:20)
est_2_10_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=14.5, β=0.25), 10., 100, 10000)), 1:20)
end

# ╔═╡ e9f6f6c5-875d-4e3a-915e-67294abacc52
#dz3, m=[0.01,0.1,1,10], u=0.01
begin
est_3_001_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=15.5, β=0.25), 0.01, 100, 10000)), 1:20)
est_3_01_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=15.5, β=0.25), 0.1, 100, 10000)), 1:20)
est_3_1_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=15.5, β=0.25), 1., 100, 10000)), 1:20)
est_3_10_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=15.5, β=0.25), 10., 100, 10000)), 1:20)
end

# ╔═╡ cab261fd-e335-4127-8872-642729aa17d3
est_3_10_001

# ╔═╡ 9eedb414-9aa1-4c64-bee6-553033e90b03
#dz4, m=[0.01,0.1,1,10], u=0.01
begin
est_4_001_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=16.5, β=0.25), 0.01, 100, 10000)), 1:20)
est_4_01_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=16.5, β=0.25), 0.1, 100, 10000)), 1:20)
est_4_1_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=16.5, β=0.25), 1., 100, 10000)), 1:20)
est_4_10_001 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UGs, μ=0., θ=16.5, β=0.25), 10., 100, 10000)), 1:20)
end

# ╔═╡ 78f588b5-8e37-4e8d-834a-0d96310717dd
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
	
#	p5_001=sum([est_5_001_001[i].p4[end] for i in 1:20] .> 100)/sum([est_5_001_001[i].pop[end] for i in 1:20] .> 100)
#	p5_01=sum([est_5_01_001[i].p4[end] for i in 1:20] .> 100)/sum([est_5_01_001[i].pop[end] for i in 1:20] .> 100)
#	p5_1=sum([est_5_1_001[i].p4[end] for i in 1:20] .> 100)/sum([est_5_1_001[i].pop[end] for i in 1:20] .> 100)
#	p5_10=sum([est_5_10_001[i].p4[end] for i in 1:20] .> 100)/sum([est_5_10_001[i].pop[end] for i in 1:20] .> 100)
	
#	p6_001=sum([est_6_001_001[i].p4[end] for i in 1:20] .> 100)/sum([est_6_001_001[i].pop[end] for i in 1:20] .> 100)
#	p6_01=sum([est_6_01_001[i].p4[end] for i in 1:20] .> 100)/sum([est_6_01_001[i].pop[end] for i in 1:20] .> 100)
#	p6_1=sum([est_6_1_001[i].p4[end] for i in 1:20] .> 100)/sum([est_6_1_001[i].pop[end] for i in 1:20] .> 100)
#	p6_10=sum([est_6_10_001[i].p4[end] for i in 1:20] .> 100)/sum([est_6_10_001[i].pop[end] for i in 1:20] .> 100)
end

# ╔═╡ b84d488c-f127-4ca3-88d2-98d2a47c68b1
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
	
#	sep5_001 = sqrt(p5_001*(1-p5_001))/20
#	sep5_01 = sqrt(p5_01*(1-p5_01))/20
#	sep5_1 = sqrt(p5_1*(1-p5_1))/20
#	sep5_10 = sqrt(p5_10*(1-p5_10))/20
	
#	sep6_001 = sqrt(p6_001*(1-p6_001))/20
#	sep6_01 = sqrt(p6_01*(1-p6_01))/20
#	sep6_1 = sqrt(p6_1*(1-p6_1))/20
#	sep6_10 = sqrt(p6_10*(1-p6_10))/20
end

# ╔═╡ 86f25465-e183-44a7-a160-596cf062ab8b
begin
estab2_01 = plot([1,2,3,4],[p1_001,p1_01,p1_1,p1_10], label="△z=1", marker = ([:hex :d]), color=:red, title="u=0.05", yerror=[sep1_001, sep1_01, sep1_1, sep1_10])
plot!([1,2,3,4],[p2_001,p2_01,p2_1,p2_10], label="△z=2", marker = ([:hex :d]), color=:orange, yerror=[sep2_001, sep2_01, sep2_1, sep2_10])
plot!([1,2,3,4],[p3_001,p3_01,p3_1,p3_10], label="△z=3", marker = ([:hex :d]), color=:yellow, yerror=[sep3_001, sep3_01, sep3_1, sep3_10])
plot!([1,2,3,4],[p4_001,p4_01,p4_1,p4_10], label="△z=4", marker = ([:hex :d]), color=:green, yerror=[sep4_001, sep4_01, sep4_1, sep4_10])
#plot!([1,2,3,4],[p5_001,p5_01,p5_1,p5_10], label="△z=5", marker = ([:hex :d]), color=:blue, ribbon=[sep5_001, sep5_01, sep5_1, sep5_10])
#plot!([1,2,3,4],[p6_001,p6_01,p6_1,p6_10], label="△z=6", marker = ([:hex :d]), color=:purple, ribbon=[sep6_001, sep6_01, sep6_1, sep6_10])
xlabel!("Migration rate")
ylabel!("Probability of tetraploid establishment")
#savefig(estab_01, "Estab_01")
end

# ╔═╡ 0543b274-d2df-4e98-a359-f27e7a8321bf
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
	
#st5_001=mean([est_5_001_001[i].ngen for i in 1:20])
#st5_01=mean([est_5_01_001[i].ngen for i in 1:20])
#st5_1=mean([est_5_1_001[i].ngen for i in 1:20])
#st5_10=mean([est_5_10_001[i].ngen for i in 1:20])
#se5_001=sem([est_5_001_001[i].ngen for i in 1:20])
#se5_01=sem([est_5_01_001[i].ngen for i in 1:20])
#se5_1=sem([est_5_1_001[i].ngen for i in 1:20])
#se5_10=sem([est_5_10_001[i].ngen for i in 1:20])
	
#st6_001=mean([est_6_001_001[i].ngen for i in 1:20])
#st6_01=mean([est_6_01_001[i].ngen for i in 1:20])
#st6_1=mean([est_6_1_001[i].ngen for i in 1:20])
#st6_10=mean([est_6_10_001[i].ngen for i in 1:20])
#se6_001=sem([est_6_001_001[i].ngen for i in 1:20])
#se6_01=sem([est_6_01_001[i].ngen for i in 1:20])
#se6_1=sem([est_6_1_001[i].ngen for i in 1:20])
#se6_10=sem([est_6_10_001[i].ngen for i in 1:20])
end

# ╔═╡ 82adacc5-8539-44e6-a688-98fd4735e691
begin
estab1_01 = plot([1,2,3,4],[st1_001, st1_01, st1_1, st1_10], label="△z=1", marker = ([:hex :d]), color=:purple, title="u=0.05", yerror=[se1_001, se1_01, se1_1, se1_10])
plot!([1,2,3,4],[st2_001, st2_01, st2_1, st2_10], label="△z=2", marker = ([:hex :d]), color=:blue, yerror=[se2_001, se2_01, se2_1, se2_10])
plot!([1,2,3,4],[st3_001, st3_01, st3_1, st3_10], label="△z=3", marker = ([:hex :d]), color=:green, yerror=[se3_001, se3_01, se3_1, se3_10])
plot!([1,2,3,4],[st4_001, st4_01, st4_1, st4_10], label="△z=4", marker = ([:hex :d]), color=:yellow, yerror=[se4_001, se4_01, se4_1, se4_10])
#plot!([1,2,3,4],[st5_001, st5_01, st5_1, st5_10], label="△z=5", marker = ([:hex :d]), color=:orange, ribbon=[se5_001, se5_01, se5_1, se5_10])
#plot!([1,2,3,4],[st6_001, st6_01, st6_1, st6_10], label="△z=6", marker = ([:hex :d]), color=:red, ribbon=[se6_001, se6_01, se6_1, se6_10])
xlabel!("Migration rate")
ylabel!("Time to establish")
#savefig(p2nsimUG, "Estab_2n_UG020020")
end

# ╔═╡ fe2ebefc-84ee-49ff-9676-57fc4e77303b
begin
	c = plot(estab1_01, estab2_01)
	#savefig(c, "Estab_005")
end

# ╔═╡ 907d41f5-1d9e-49ef-be56-7313433525f8
plot([1,2,3,4],[st1_001, st1_01, st1_1, st1_10], label="△z=1", marker = ([:hex :d]), color=:green, title="u=0.01", ribbon=[se1_001, se1_01, se1_1, se1_10])

# ╔═╡ Cell order:
# ╟─9aa821e0-a906-11eb-10e2-89580ee28b2c
# ╠═fb6945a6-5c05-474d-b8f8-2749fd304a30
# ╠═e9d83c94-d3c4-41ca-b4d8-30ee1eda413d
# ╠═5cfb3522-997e-42c3-8156-26fe6bf615c0
# ╠═41ffb48d-ffdf-4b24-bbf6-4c3194295d71
# ╠═e9f6f6c5-875d-4e3a-915e-67294abacc52
# ╠═9eedb414-9aa1-4c64-bee6-553033e90b03
# ╠═3609e21d-da26-4cb6-93fe-d3632d19d3d7
# ╠═95c2aa35-fb67-4dea-9cf0-d7d8d8b4d31c
# ╠═78f588b5-8e37-4e8d-834a-0d96310717dd
# ╠═b84d488c-f127-4ca3-88d2-98d2a47c68b1
# ╠═cab261fd-e335-4127-8872-642729aa17d3
# ╠═86f25465-e183-44a7-a160-596cf062ab8b
# ╠═0543b274-d2df-4e98-a359-f27e7a8321bf
# ╠═82adacc5-8539-44e6-a688-98fd4735e691
# ╠═5e04a741-c863-4689-95f0-e55370e53d7d
# ╠═fe2ebefc-84ee-49ff-9676-57fc4e77303b
# ╠═907d41f5-1d9e-49ef-be56-7313433525f8
# ╟─b690a6ed-5e38-40b4-ac80-e2f9ed3882c4
# ╠═f6f9d8d9-e253-4431-9690-802511aa869f
# ╠═d4c6dd9e-b993-49cc-910b-7d408c5fe500
# ╟─adbeed31-db8c-46ce-9c84-ee48d2173c9c
