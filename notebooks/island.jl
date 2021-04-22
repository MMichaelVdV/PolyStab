### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 1f791df0-8b5e-11eb-3bc5-e7302263af09
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, StatsPlots

# ╔═╡ 325d63e0-8b5e-11eb-0d55-f149adc308e1
using PolyStab

# ╔═╡ 423cffe0-8b5f-11eb-1df4-ff7eed10289a
using PolyStab: Agent, randagent_p, MixedPloidyDeme, IslandDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, random_mating, evolving_neutraldeme, evolving_islanddeme

# ╔═╡ 44771223-a589-4de3-9fea-321b70138d80
using PolyStab: malthusian_fitness, trait_mean, ploidy_freq, f_trait_agents, mating_PnB_x

# ╔═╡ a7bc4e02-8b5c-11eb-0e36-ef47de093003
md""" ### Island model
"""

# ╔═╡ ba38834e-8b5c-11eb-246c-99429ab2ae27
md""" ###### cfr. Establishment in a new habitat by polygenic adaptation (Barton and Etheridge, 2017)"""

# ╔═╡ baffcaf2-8b5c-11eb-1479-e110e25be9d5
md""" Representation: large abstract mainland population (assumes: infinite population size in HWLE), migration from mainland to (different smaller) island(s) with finite population size (can be respresented as demes).

Some quick thoughts:
- mainland -> island: establishment of polyploid in a diploid island population
- migration between islands (demes) or only from mainland to island?
- starting population mainland (Gv), island
- colonization event (selfing versus non-selfing, i.e. with nonselfing you need at least two individuals to succesfully colonize an island). 
- stochasticity < dispersal (phenotype - island optimum mismatch), i.e. do all islands have the same optimum as mainland (agents will already be adapted) or is there variation in island optima (different islands with different parameters versus all islands same parameters and different simulation runs).
"""

# ╔═╡ 164cbf80-8b5d-11eb-129c-2fd4863912df
md""" ###### Some simulations"""

# ╔═╡ 3712a3fe-8b5d-11eb-2560-5940bd21ea0c
md""" First we can look a the case with a single migrant migrating to a single island. Simulations start with an idealized infinite mainland population in HWLE. There is migration from this population to the island with migration rate M. Selfing is allowed so a single migrant can establish a new population on an island. The source is assumed to be poorly adapted to conditions on the island with selection gradient β."""

# ╔═╡ d97a6a70-8b5d-11eb-0ed6-056d834ce9bb
md""" The island can be modelled as a single deme:"""

# ╔═╡ 2625fd72-8b5f-11eb-1e69-adf12cbd84c9
island = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.)

# ╔═╡ c4067530-8fcc-11eb-322e-bdeba99f75a4
islanddeme = IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], β=.2)

# ╔═╡ 6c24a8d9-7ecd-471a-b3d8-7b7ab6ef32d2
begin
n = 50. #number of loci
α = 0.5 #allelic effect size
d = 1
E_mean2 = n*α*d*0.5
E_var2 = 0.5 * n * (α^2 / 4)
distro = map(x->trait(randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1)[1]),1:500)
	
plot(Normal(E_mean2,sqrt(E_var2)),color=:black, label="Expected 2N")
plot!(Normal(mean(distro), std(distro)), label="Observed 2N")
xlabel!("\$phenotype\$")
ylabel!("\$\$")
end

# ╔═╡ 8272fcde-8b5f-11eb-1c13-b591229f902f
migrant = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1)[1]

# ╔═╡ 93754ff5-0428-4188-bf46-49137ee1fc4f
trait(migrant)

# ╔═╡ faa45f90-8b6b-11eb-2826-776161a15a58
migrant.loci

# ╔═╡ ce694600-8bd6-11eb-0d50-c7bd17294608
allelefreqs_p(migrant)

# ╔═╡ 05a7f6f0-8bd9-11eb-0892-f3cde611f9ee
mean(allelefreqs_p(migrant))

# ╔═╡ 143183d0-8b61-11eb-1986-310271688189
push!(island.agents,migrant)

# ╔═╡ fc0f65e0-8fcc-11eb-21ab-dd9a09e00355
push!(islanddeme.agents,migrant)

# ╔═╡ 4d1ba2c0-8bcf-11eb-2b00-adaff59d3e9e
length(island.agents)

# ╔═╡ 7b190320-8bcf-11eb-15b6-2b3a871bc3bc
allelefreqs_p(island)

# ╔═╡ 93164980-8bcd-11eb-02a1-496109bc20c3
2 .*heterozygosities_p(island)

# ╔═╡ 91834c70-8bd8-11eb-0ef6-e79beaa1849b
mean(2 .*heterozygosities_p(island))

# ╔═╡ 52c76560-8b61-11eb-1edd-e71d92ca9a1e
island

# ╔═╡ 08db9cd0-8fcd-11eb-1e0a-fd10a73b4dd0
islanddeme

# ╔═╡ 58605950-8b61-11eb-16d9-0ff9be46d489
island_p = evolving_selectiondeme(island,20)

# ╔═╡ 145aba50-8fcd-11eb-081d-4ffe44449e3f
island_p2 = evolving_islanddeme(islanddeme,20)

# ╔═╡ 6df662b0-8b6a-11eb-0e74-611b7747413e
begin
	pf2_p1 = island_p.p2
	pf3_p1 = island_p.p3
	pf4_p1 = island_p.p4
	p1 = plot(pf2_p1, grid=false, color=:blue, label="diploids", legend=:bottomright)
	#plot!(pf3_p1, grid=false, color=:green, label="triploids")
	plot!(pf4_p1, grid=false, color=:red, label="tetraploids")
	hline!([island.K],label ="K",colour = "black",linestyle=:dash)
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ dc3ad710-8b6a-11eb-3e55-7ba4760f8132
begin
	traitmean_ploidy = map(mean, island_p.tm)
	p2 = plot(traitmean_ploidy, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright) #title="Stabilizing selection"
	for (i,t) in enumerate(island_p.fta)
	scatter!([i for x in 1:length(t)],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([island.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ 84585d6e-8c09-11eb-3c5c-af38eb556b15
var_add(a::Agent,α) = ploidy(a)*α^2*sum(heterozygosities_p(a))

# ╔═╡ 7473dcf0-8fcd-11eb-16e8-a5371cf5105f
begin
	pf2_p2 = island_p2.p2
	pf3_p2 = island_p2.p3
	pf4_p2 = island_p2.p4
	p21 = plot(pf2_p2, grid=false, color=:blue, label="diploids", legend=:bottomright)
	#plot!(pf3_p1, grid=false, color=:green, label="triploids")
	plot!(pf4_p2, grid=false, color=:red, label="tetraploids")
	hline!([island.K],label ="K",colour = "black",linestyle=:dash)
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ 1bbe42d0-8fd2-11eb-0d53-8be7e42d7722
begin
	traitmean_ploidy2 = map(mean, island_p2.tm)
	p22 = plot(traitmean_ploidy2, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright) #title="Stabilizing selection"
	for (i,t) in enumerate(island_p2.fta)
	scatter!([i for x in 1:length(t)],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([island.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ 1cd88c70-8fd2-11eb-30bc-a3b3c916b97d
begin
	Hₒ_ploidyhet2 = map(mean, 2 .*island_p2.het)
	expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	p23 = plot(Hₒ_ploidyhet2, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH")
	plot!(1:island_p2.ngen, 
		t->expected_heterozygosity(Hₒ_ploidyhet2[1], t, 50),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	V_A2 = map(sum, 2*(0.1)^2 .*island_p2.het)
	plot!(V_A2, grid=false, color=:red, label="\$V_A(t)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ ee422a80-8b6a-11eb-13ce-a7f8ab3de0ef
begin
	Hₒ_ploidyhet = map(mean, 2 .*island_p.het)
	p3 = plot(Hₒ_ploidyhet, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH")
	plot!(1:island_p.ngen, 
		t->expected_heterozygosity(Hₒ_ploidyhet[1], t, 50),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	V_A = map(sum, 2*(0.1)^2 .*island_p.het)
	plot!(V_A, grid=false, color=:red, label="\$V_A(t)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 5a0f4c50-8c26-11eb-0f10-dfd79a2ff7c5
begin
	plot(p1,p2,p3)
	#savefig("single_migrant.pdf")
end

# ╔═╡ 4be009d0-8fd2-11eb-2338-15f33ab00a20
plot(p21,p22,p23)

# ╔═╡ 4ce3fdbb-82f2-46da-a226-485306b0d91b


# ╔═╡ 602cebcc-46bf-48fe-a3d3-163318796b98
md""" ### Simulations"""

# ╔═╡ 1182e78a-f422-4e14-b6f4-e6f2c05c3da4
island_s2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.)

# ╔═╡ c3814ccf-976e-4c9c-99dc-ddab6b4945b8
#migrant_s2 = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1)[1]

# ╔═╡ 9db991e8-49af-44be-abdb-646d5409d879
#trait(migrant_s2)

# ╔═╡ 6b38a600-768f-4abc-9e3f-63a010b15461
#push!(island_s2.agents,migrant_s2)

# ╔═╡ 8841a695-1110-487e-bd8a-720ef80a64c3
s2 = map(x->(evolving_selectiondeme(MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.),20)).pop[end], 1:1000) 

# ╔═╡ b3f9e34a-2a82-4223-ac2a-6c9e8442b0ab
begin
	histogram(s2, bins = 25, fillalpha = 0.4, title="Estab of diploid migrant")
	xlabel!("\$popsize\$")
	ylabel!("\$n\$")
end

# ╔═╡ 1265d278-030c-48f9-8f2b-b73375b64823
length(s2[s2 .> 0])/1000

# ╔═╡ 06a553ee-662d-4bbf-8a5f-365912b14175
island_s4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.)

# ╔═╡ d7dae018-84ac-493d-96a3-d400073b29a2
#migrant_s4 = randagent_p(0.5, 0.5, 50, [0., 0., 0., 1.],1)[1]

# ╔═╡ 2d4c5f74-07a5-4082-b75b-93a772376173
#trait(migrant_s4)

# ╔═╡ 65cbc68c-7176-4bc1-8a2a-20f245311918
#push!(island_s4.agents,migrant_s4)

# ╔═╡ 0d4614ce-3d3b-4b57-ba53-9b42a11c711b
s4 = map(x->(evolving_selectiondeme(MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 0., 0., 1.],1), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.),20)).pop[end], 1:1000) 

# ╔═╡ 091ec75c-9d54-4b3d-bc45-dbb748ca447a
begin
	histogram(s4, bins = 25, fillalpha = 0.4, title="Estab of tetraploid migrant")
	xlabel!("\$popsize\$")
	ylabel!("\$n\$")
end

# ╔═╡ 558da9c6-e3f1-481e-86c2-8b825ae1f75f
length(s4[s4 .> 0])/1000

# ╔═╡ 71012460-b6af-4282-a78d-692d0a2a9577
md""" IMPORTANT REMARK: The starting phenotypic variance of the island simulated this way for 4N is only half of that of 2N, i.e. there will be less migrants with extreme phenotypes. """

# ╔═╡ 318e6c52-4a13-40af-bc12-af9b1185e6f7
md""" ### Continuous migration """

# ╔═╡ 95ee86cd-14df-42da-b7b6-e094551c298e
function evolving_islandwmigration(d::MixedPloidyDeme, ngen; 
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
		mig = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1)[1]
		push!(d.agents,mig)
		d = mating_PnB_x(d)
		#d = mutate(d) 
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

# ╔═╡ c60e7da5-3294-4ca0-908c-d3941b559c62
island_m2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.)

# ╔═╡ 5645c34b-5049-49e7-ab47-c0cb33e3e6a8
sim_mig = evolving_islandwmigration(island_m2, 20)

# ╔═╡ f1595c57-e039-41aa-b8e9-9997cdc19641
begin
	pf2_p2m = sim_mig.p2
	pf3_p2m = sim_mig.p3
	pf4_p2m = sim_mig.p4
	p21m = plot(pf2_p2m, grid=false, color=:blue, label="diploids", legend=:bottomright)
	#plot!(pf3_p1, grid=false, color=:green, label="triploids")
	plot!(pf4_p2m, grid=false, color=:red, label="tetraploids")
	hline!([island.K],label ="K",colour = "black",linestyle=:dash)
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ ce63121d-7c50-49d2-a60c-1771bcd01ae9
begin
	traitmean_mig = map(mean, sim_mig.tm)
	p2m = plot(traitmean_mig, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright) #title="Stabilizing selection"
	for (i,t) in enumerate(sim_mig.fta)
	scatter!([i for x in 1:length(t)],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([island.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ f51ad793-1dc2-4515-8b7c-627904bc016e
begin
	Hₒ_ploidyhet2m = map(mean, 2 .*sim_mig.het)
	#expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	p23m = plot(Hₒ_ploidyhet2m, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH")
	plot!(1:sim_mig.ngen, 
		t->expected_heterozygosity(Hₒ_ploidyhet2m[1], t, 50),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	V_A2m = map(sum, 2*(0.1)^2 .*sim_mig.het)
	plot!(V_A2, grid=false, color=:red, label="\$V_A(t)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ Cell order:
# ╟─a7bc4e02-8b5c-11eb-0e36-ef47de093003
# ╟─ba38834e-8b5c-11eb-246c-99429ab2ae27
# ╟─baffcaf2-8b5c-11eb-1479-e110e25be9d5
# ╟─164cbf80-8b5d-11eb-129c-2fd4863912df
# ╠═1f791df0-8b5e-11eb-3bc5-e7302263af09
# ╠═325d63e0-8b5e-11eb-0d55-f149adc308e1
# ╠═423cffe0-8b5f-11eb-1df4-ff7eed10289a
# ╟─3712a3fe-8b5d-11eb-2560-5940bd21ea0c
# ╟─d97a6a70-8b5d-11eb-0ed6-056d834ce9bb
# ╠═2625fd72-8b5f-11eb-1e69-adf12cbd84c9
# ╠═c4067530-8fcc-11eb-322e-bdeba99f75a4
# ╠═6c24a8d9-7ecd-471a-b3d8-7b7ab6ef32d2
# ╠═8272fcde-8b5f-11eb-1c13-b591229f902f
# ╠═93754ff5-0428-4188-bf46-49137ee1fc4f
# ╠═faa45f90-8b6b-11eb-2826-776161a15a58
# ╠═ce694600-8bd6-11eb-0d50-c7bd17294608
# ╠═05a7f6f0-8bd9-11eb-0892-f3cde611f9ee
# ╠═143183d0-8b61-11eb-1986-310271688189
# ╠═fc0f65e0-8fcc-11eb-21ab-dd9a09e00355
# ╠═4d1ba2c0-8bcf-11eb-2b00-adaff59d3e9e
# ╠═7b190320-8bcf-11eb-15b6-2b3a871bc3bc
# ╠═93164980-8bcd-11eb-02a1-496109bc20c3
# ╠═91834c70-8bd8-11eb-0ef6-e79beaa1849b
# ╠═52c76560-8b61-11eb-1edd-e71d92ca9a1e
# ╠═08db9cd0-8fcd-11eb-1e0a-fd10a73b4dd0
# ╠═58605950-8b61-11eb-16d9-0ff9be46d489
# ╠═145aba50-8fcd-11eb-081d-4ffe44449e3f
# ╠═6df662b0-8b6a-11eb-0e74-611b7747413e
# ╠═dc3ad710-8b6a-11eb-3e55-7ba4760f8132
# ╠═ee422a80-8b6a-11eb-13ce-a7f8ab3de0ef
# ╠═84585d6e-8c09-11eb-3c5c-af38eb556b15
# ╠═7473dcf0-8fcd-11eb-16e8-a5371cf5105f
# ╠═1bbe42d0-8fd2-11eb-0d53-8be7e42d7722
# ╠═1cd88c70-8fd2-11eb-30bc-a3b3c916b97d
# ╠═5a0f4c50-8c26-11eb-0f10-dfd79a2ff7c5
# ╠═4be009d0-8fd2-11eb-2338-15f33ab00a20
# ╠═4ce3fdbb-82f2-46da-a226-485306b0d91b
# ╟─602cebcc-46bf-48fe-a3d3-163318796b98
# ╠═1182e78a-f422-4e14-b6f4-e6f2c05c3da4
# ╠═c3814ccf-976e-4c9c-99dc-ddab6b4945b8
# ╠═9db991e8-49af-44be-abdb-646d5409d879
# ╠═6b38a600-768f-4abc-9e3f-63a010b15461
# ╠═8841a695-1110-487e-bd8a-720ef80a64c3
# ╠═b3f9e34a-2a82-4223-ac2a-6c9e8442b0ab
# ╠═1265d278-030c-48f9-8f2b-b73375b64823
# ╠═06a553ee-662d-4bbf-8a5f-365912b14175
# ╠═d7dae018-84ac-493d-96a3-d400073b29a2
# ╠═2d4c5f74-07a5-4082-b75b-93a772376173
# ╠═65cbc68c-7176-4bc1-8a2a-20f245311918
# ╠═0d4614ce-3d3b-4b57-ba53-9b42a11c711b
# ╠═091ec75c-9d54-4b3d-bc45-dbb748ca447a
# ╠═558da9c6-e3f1-481e-86c2-8b825ae1f75f
# ╟─71012460-b6af-4282-a78d-692d0a2a9577
# ╟─318e6c52-4a13-40af-bc12-af9b1185e6f7
# ╠═44771223-a589-4de3-9fea-321b70138d80
# ╠═95ee86cd-14df-42da-b7b6-e094551c298e
# ╠═c60e7da5-3294-4ca0-908c-d3941b559c62
# ╠═5645c34b-5049-49e7-ab47-c0cb33e3e6a8
# ╠═f1595c57-e039-41aa-b8e9-9997cdc19641
# ╠═ce63121d-7c50-49d2-a60c-1771bcd01ae9
# ╠═f51ad793-1dc2-4515-8b7c-627904bc016e
