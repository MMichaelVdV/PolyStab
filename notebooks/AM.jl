### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 88d20712-ab37-49d7-9baf-5065697aac64
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, DataFrames, StatsPlots

# ╔═╡ dcb3fcf0-8769-4335-8bb8-6801f68e7b57
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, evolving_neutraldeme, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy, trait_mean, randagent, evolving_haploiddeme, mate_p, evolving_selectiondemeh, malthusian_fitness, number_of_offspring, ismock, ploidy_freq, f_trait_agents, mutate

# ╔═╡ e8c1c400-d78f-11eb-2a0b-2914730bafcb
md"""### Assortative mating"""

# ╔═╡ b4eabf0b-c15d-4e02-8bf0-cf7afe3c6a45
md""" H: Assortative mating can help to overcome minority cytotype exclusion. Effects on inbreeding depression, effects in finite population size (drift?), ...""" 

# ╔═╡ 5b99efb0-ab9b-42e7-bf84-c1fb1787669a
d_p1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [0., 1., 0., 0.], 150), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 0.1)

# ╔═╡ 2dcaedf3-c1f6-4170-a848-bf1c9e3d89ce
ploidy.(d_p1.agents)

# ╔═╡ c63bb0d3-b881-4d1e-9766-1f41991aceba
"""
	mating_PnB_s(d::MixedPloidyDeme{A})

Mating in a mixed ploidy deme with unreduced gamete formation and partner
choice weighted by malthusian fitness.
"""
function mating_PnB_a(d::AbstractDeme{A}, a::Float64) where A
	new_agents = A[]
	fitnesses = exp.(malthusian_fitness(d))
	pl = ploidy.(d.agents)
	for i=1:length(d)
		B1 = d.agents[i]
		if rand() < a 
			pls = (pl .*0) .+ (pl.==ploidy(B1))
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

# ╔═╡ b6d363f2-5954-410f-b750-36b725f99a66
"""
	mating_PnB_s(d::MixedPloidyDeme{A})

Mating in a mixed ploidy deme with unreduced gamete formation and partner
choice weighted by malthusian fitness.
"""
function mating_PnB_s(d::AbstractDeme{A}, s::Float64) where A
	new_agents = A[]
	fitnesses = exp.(malthusian_fitness(d))
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

# ╔═╡ d21328d6-6496-45d1-954c-e0921c013eb6
trait.(d_p1.agents)

# ╔═╡ ccf77978-18ed-41db-9bb8-eeee2e4013f8
ploidy.(d_p1.agents) .== ploidy(d_p1.agents[1])

# ╔═╡ 04bbf67f-6c8a-4e07-816b-3c4c724b8241
fitnesses = exp.(malthusian_fitness(d_p1))

# ╔═╡ 4f56ee9b-d8f8-4c1a-b7fd-5036509e5e8f
weights(fitnesses)

# ╔═╡ 14e77135-bd00-4818-acbc-75d2c48cf980
sum(weights(fitnesses))

# ╔═╡ f1f4216d-cd11-41d6-a4cd-36a292958167
mating_PnB_s(d_p1, 0.1)

# ╔═╡ d3eac2a7-ff67-418d-b6f1-0fb4e32f7868
"""
	evolving_selectiondeme(d::AbstractDeme, ngen)

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.
"""
function evolving_selectiondeme_s(d::MixedPloidyDeme, s, ngen; 
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
		d = mating_PnB_s(d, s)
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

# ╔═╡ 7af043fa-36ed-418b-a523-eceee16f9d55
begin
	sel_p2_00 = evolving_selectiondeme_s(d_p1, 0., 1000)
	sel_p2_02 = evolving_selectiondeme_s(d_p1, 0.2, 1000)
	sel_p2_04 = evolving_selectiondeme_s(d_p1, 0.4, 1000)
	sel_p2_06 = evolving_selectiondeme_s(d_p1, 0.6, 1000)
	sel_p2_08 = evolving_selectiondeme_s(d_p1, 0.8, 1000)
	sel_p2_10 = evolving_selectiondeme_s(d_p1, 1., 1000)
end

# ╔═╡ 158f286c-7c34-48de-a56d-d7ec177ad64c
begin
	Hₒ_sel_p2_00 = map(mean, sel_p2_00.het)
	Hₒ_sel_p2_02 = map(mean, sel_p2_02.het)
	Hₒ_sel_p2_04 = map(mean, sel_p2_04.het)
	Hₒ_sel_p2_06 = map(mean, sel_p2_06.het)
	Hₒ_sel_p2_08 = map(mean, sel_p2_08.het)
	Hₒ_sel_p2_10 = map(mean, sel_p2_10.het)


	
	plot(Hₒ_sel_p2_00, grid=false, color=:black, label="\$H_op1(t)\$", title="LOH")
	plot!(Hₒ_sel_p2_02, grid=false, color=:blue, label="\$H_op2(t)\$")
	plot!(Hₒ_sel_p2_04, grid=false, color=:blue, label="\$H_op2(t)\$")
	plot!(Hₒ_sel_p2_06, grid=false, color=:blue, label="\$H_op2(t)\$")
	plot!(Hₒ_sel_p2_08, grid=false, color=:blue, label="\$H_op2(t)\$")
	plot!(Hₒ_sel_p2_10, grid=false, color=:blue, label="\$H_op2(t)\$")	
	

	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 5b0f78ca-b778-4219-8d30-2a47fcf1b13b
begin
	trait_p2_00 = map(mean, sel_p2_00.tm)
	p2n_00 = plot(trait_p2_00, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid", legend=:bottomright)
	for (i,t) in enumerate(sel_p2_00.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(17,25)
end

# ╔═╡ 66d0a754-e4de-4b83-a4b3-6def51d80284
begin
	trait_p2_02 = map(mean, sel_p2_02.tm)
	p2n_02 = plot(trait_p2_02, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid", legend=:bottomright)
	for (i,t) in enumerate(sel_p2_02.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(17,25)
end

# ╔═╡ b313d866-28c9-41b7-a765-438bfcc747b6
begin
	trait_p2_04 = map(mean, sel_p2_04.tm)
	p2n_04 = plot(trait_p2_04, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid", legend=:bottomright)
	for (i,t) in enumerate(sel_p2_04.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(17,25)
end

# ╔═╡ bb562ff0-b950-4921-8daf-0052d4e7d3ff
begin
	trait_p2_06 = map(mean, sel_p2_06.tm)
	p2n_06 = plot(trait_p2_06, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid", legend=:bottomright)
	for (i,t) in enumerate(sel_p2_06.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(17,25)
end

# ╔═╡ 5fb215cb-59e9-4088-a868-d072cfbdea56
begin
	trait_p2_08 = map(mean, sel_p2_08.tm)
	p2n_08 = plot(trait_p2_08, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid", legend=:bottomright)
	for (i,t) in enumerate(sel_p2_08.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(17,25)
end

# ╔═╡ 3a9601db-1075-48cd-9eec-1e532a0e920a
begin
	trait_p2_10 = map(mean, sel_p2_10.tm)
	p2n_10 = plot(trait_p2_10, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid", legend=:bottomright)
	for (i,t) in enumerate(sel_p2_10.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(17,25)
end

# ╔═╡ 7c57e021-f347-4900-ace0-1a1c4fbe59f3
begin
	pop_p2_00 = map(mean, sel_p2_00.p2)
	plot(pop_p2_00, grid=false, color=:blue, label="\$pop_p2(t)\$")
	pop_p2_02 = map(mean, sel_p2_02.p2)
	plot!(pop_p2_02, grid=false, color=:blue, label="\$pop_p2(t)\$")
	pop_p2_04 = map(mean, sel_p2_04.p2)
	plot!(pop_p2_04, grid=false, color=:blue, label="\$pop_p2(t)\$")
	pop_p2_06 = map(mean, sel_p2_06.p2)
	plot!(pop_p2_06, grid=false, color=:blue, label="\$pop_p2(t)\$")
	pop_p2_08 = map(mean, sel_p2_08.p2)
	plot!(pop_p2_08, grid=false, color=:blue, label="\$pop_p2(t)\$")
	pop_p2_10 = map(mean, sel_p2_10.p2)
	plot!(pop_p2_10, grid=false, color=:blue, label="\$pop_p2(t)\$")
	xlabel!("\$t\$")
	#ylims!(145,155)
	ylabel!("\$pop(t)\$")
end

# ╔═╡ b733ac9a-a34b-4301-9da8-a9c4acf09b46
begin
	pop_p4_00 = map(mean, sel_p2_00.p4)
	plot(pop_p4_00, grid=false, color=:blue, label="\$pop_p2(t)\$")
	pop_p4_02 = map(mean, sel_p2_02.p4)
	plot!(pop_p4_02, grid=false, color=:blue, label="\$pop_p2(t)\$")
	pop_p4_04 = map(mean, sel_p2_04.p4)
	plot!(pop_p4_04, grid=false, color=:blue, label="\$pop_p2(t)\$")
	pop_p4_06 = map(mean, sel_p2_06.p4)
	plot!(pop_p4_06, grid=false, color=:blue, label="\$pop_p2(t)\$")
	pop_p4_08 = map(mean, sel_p2_08.p4)
	plot!(pop_p4_08, grid=false, color=:blue, label="\$pop_p2(t)\$")
	pop_p4_10 = map(mean, sel_p2_10.p4)
	plot!(pop_p4_10, grid=false, color=:blue, label="\$pop_p2(t)\$")
	xlabel!("\$t\$")
	#ylims!(145,155)
	ylabel!("\$pop(t)\$")
end

# ╔═╡ d91ad275-e595-43f7-aebc-103d8adb3787
md""" Mixed ploidy scenario """

# ╔═╡ Cell order:
# ╟─e8c1c400-d78f-11eb-2a0b-2914730bafcb
# ╟─b4eabf0b-c15d-4e02-8bf0-cf7afe3c6a45
# ╠═88d20712-ab37-49d7-9baf-5065697aac64
# ╠═dcb3fcf0-8769-4335-8bb8-6801f68e7b57
# ╠═5b99efb0-ab9b-42e7-bf84-c1fb1787669a
# ╠═2dcaedf3-c1f6-4170-a848-bf1c9e3d89ce
# ╠═c63bb0d3-b881-4d1e-9766-1f41991aceba
# ╠═b6d363f2-5954-410f-b750-36b725f99a66
# ╠═d21328d6-6496-45d1-954c-e0921c013eb6
# ╠═ccf77978-18ed-41db-9bb8-eeee2e4013f8
# ╠═04bbf67f-6c8a-4e07-816b-3c4c724b8241
# ╠═4f56ee9b-d8f8-4c1a-b7fd-5036509e5e8f
# ╠═14e77135-bd00-4818-acbc-75d2c48cf980
# ╠═f1f4216d-cd11-41d6-a4cd-36a292958167
# ╠═d3eac2a7-ff67-418d-b6f1-0fb4e32f7868
# ╠═7af043fa-36ed-418b-a523-eceee16f9d55
# ╠═158f286c-7c34-48de-a56d-d7ec177ad64c
# ╟─5b0f78ca-b778-4219-8d30-2a47fcf1b13b
# ╟─66d0a754-e4de-4b83-a4b3-6def51d80284
# ╟─b313d866-28c9-41b7-a765-438bfcc747b6
# ╟─bb562ff0-b950-4921-8daf-0052d4e7d3ff
# ╟─5fb215cb-59e9-4088-a868-d072cfbdea56
# ╟─3a9601db-1075-48cd-9eec-1e532a0e920a
# ╠═7c57e021-f347-4900-ace0-1a1c4fbe59f3
# ╠═b733ac9a-a34b-4301-9da8-a9c4acf09b46
# ╠═d91ad275-e595-43f7-aebc-103d8adb3787
