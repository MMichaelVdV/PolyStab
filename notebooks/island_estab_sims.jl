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
us = 0.01

# ╔═╡ fb6945a6-5c05-474d-b8f8-2749fd304a30
UGs = [0. 0. 0. 0. ; 1-us us 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.]

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
end

# ╔═╡ 86f25465-e183-44a7-a160-596cf062ab8b
begin
estab2_01 = plot([1,2,3,4],[p1_001,p1_01,p1_1,p1_10], label="△z=1", marker = ([:hex :d]), color=:green, title="u=0.01")
plot!([1,2,3,4],[p2_001,p2_01,p2_1,p2_10], label="△z=2", marker = ([:hex :d]), color=:blue)
plot!([1,2,3,4],[p3_001,p3_01,p3_1,p3_10], label="△z=3", marker = ([:hex :d]), color=:red)
plot!([1,2,3,4],[p4_001,p4_01,p4_1,p4_10], label="△z=4", marker = ([:hex :d]), color=:black)
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

st2_001=mean([est_2_001_001[i].ngen for i in 1:20])
st2_01=mean([est_2_01_001[i].ngen for i in 1:20])
st2_1=mean([est_2_1_001[i].ngen for i in 1:20])
st2_10=mean([est_2_10_001[i].ngen for i in 1:20])

st3_001=mean([est_3_001_001[i].ngen for i in 1:20])
st3_01=mean([est_3_01_001[i].ngen for i in 1:20])
st3_1=mean([est_3_1_001[i].ngen for i in 1:20])
st3_10=mean([est_3_10_001[i].ngen for i in 1:20])

st4_001=mean([est_4_001_001[i].ngen for i in 1:20])
st4_01=mean([est_4_01_001[i].ngen for i in 1:20])
st4_1=mean([est_4_1_001[i].ngen for i in 1:20])
st4_10=mean([est_4_10_001[i].ngen for i in 1:20])
end

# ╔═╡ 82adacc5-8539-44e6-a688-98fd4735e691
begin
estab1_01 = plot([1,2,3,4],[st1_001, st1_01, st1_1, st1_10], label="△z=1", marker = ([:hex :d]), color=:green, title="u=0.01")
plot!([1,2,3,4],[st2_001, st2_01, st2_1, st2_10], label="△z=2", marker = ([:hex :d]), color=:blue)
plot!([1,2,3,4],[st3_001, st3_01, st3_1, st3_10], label="△z=3", marker = ([:hex :d]), color=:red)
plot!([1,2,3,4],[st4_001, st4_01, st4_1, st4_10], label="△z=4", marker = ([:hex :d]), color=:black)
xlabel!("Migration rate")
ylabel!("Time to establish")
#savefig(p2nsimUG, "Estab_2n_UG020020")
end

# ╔═╡ fe2ebefc-84ee-49ff-9676-57fc4e77303b
begin
	c = plot(estab1_01, estab2_01)
	savefig(c, "Estab_001")
end

# ╔═╡ Cell order:
# ╟─9aa821e0-a906-11eb-10e2-89580ee28b2c
# ╠═fb6945a6-5c05-474d-b8f8-2749fd304a30
# ╠═e9d83c94-d3c4-41ca-b4d8-30ee1eda413d
# ╠═5cfb3522-997e-42c3-8156-26fe6bf615c0
# ╠═41ffb48d-ffdf-4b24-bbf6-4c3194295d71
# ╠═e9f6f6c5-875d-4e3a-915e-67294abacc52
# ╠═9eedb414-9aa1-4c64-bee6-553033e90b03
# ╠═78f588b5-8e37-4e8d-834a-0d96310717dd
# ╠═cab261fd-e335-4127-8872-642729aa17d3
# ╠═86f25465-e183-44a7-a160-596cf062ab8b
# ╠═0543b274-d2df-4e98-a359-f27e7a8321bf
# ╠═82adacc5-8539-44e6-a688-98fd4735e691
# ╠═fe2ebefc-84ee-49ff-9676-57fc4e77303b
# ╟─b690a6ed-5e38-40b4-ac80-e2f9ed3882c4
# ╠═f6f9d8d9-e253-4431-9690-802511aa869f
# ╠═d4c6dd9e-b993-49cc-910b-7d408c5fe500
# ╟─adbeed31-db8c-46ce-9c84-ee48d2173c9c
