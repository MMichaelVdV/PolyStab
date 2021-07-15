### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 1f791df0-8b5e-11eb-3bc5-e7302263af09
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, StatsPlots

# ╔═╡ 325d63e0-8b5e-11eb-0d55-f149adc308e1
using PolyStab

# ╔═╡ 423cffe0-8b5f-11eb-1df4-ff7eed10289a
using PolyStab: Agent, randagent_p, MixedPloidyDeme, IslandDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, random_mating, evolving_neutraldeme, evolving_islanddeme

# ╔═╡ 44771223-a589-4de3-9fea-321b70138d80
using PolyStab: malthusian_fitness, trait_mean, ploidy_freq, f_trait_agents, mating_PnB_x, directional_selection, mating_PnB, mutate

# ╔═╡ c0a08e7f-dca4-4107-b2b3-d18d7dae77a8
using PolyStab: popsize

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
island = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.)

# ╔═╡ c4067530-8fcc-11eb-322e-bdeba99f75a4
islanddeme = IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], β=.25, θ=13.)

# ╔═╡ bfc7fa55-15a7-4690-9c33-7cfffedc975c
islanddeme_UG = IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.9 0.1 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], β=.25, θ=13.)

# ╔═╡ 6c24a8d9-7ecd-471a-b3d8-7b7ab6ef32d2
begin
n = 250. #number of loci
α = 0.1 #allelic effect size
d = 1
E_mean2 = n*α*d*0.5
E_var2 = 0.5 * n * (α^2 / 4)
distro = map(x->trait(randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],1)[1]),1:500)
	
plot(Normal(E_mean2,sqrt(E_var2)),color=:black, label="Expected 2N")
plot!(Normal(mean(distro), std(distro)), label="Observed 2N")
xlabel!("\$phenotype\$")
ylabel!("\$\$")
end

# ╔═╡ 09247474-d2d7-4a19-b3b6-dccf77fcd7f0
sort(distro)

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
island_p = evolving_selectiondeme(island,30)

# ╔═╡ 145aba50-8fcd-11eb-081d-4ffe44449e3f
island_p2 = evolving_islanddeme(islanddeme,30)

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
 md""" ### Black-hole sink dynamics"""

# ╔═╡ e87406ec-2deb-4fda-9e49-5f95676b32c3
md""" A black-hole sink a metaphor for a deme with one-way immigration from a source population and without emigration."""

# ╔═╡ 602cebcc-46bf-48fe-a3d3-163318796b98
 md""" ### Simulations with a single migrant"""

# ╔═╡ 1182e78a-f422-4e14-b6f4-e6f2c05c3da4
island_s2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.)

# ╔═╡ 824018c1-6f6b-40ff-8b77-aebaccc20afc
island_s2i = IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15., β=0.25)

# ╔═╡ 8841a695-1110-487e-bd8a-720ef80a64c3
s2 = map(x->(evolving_selectiondeme(MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.),20)).pop[end], 1:1000) 

# ╔═╡ 2da9fafc-64c4-4c4c-a0c8-2e248ae84441
s2i = map(x->(evolving_islanddeme(IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15., β=0.25),20)).pop[end], 1:1000) 

# ╔═╡ ef4e8ee0-828f-478c-ba12-1a555d620019
sort(s2i)

# ╔═╡ b3f9e34a-2a82-4223-ac2a-6c9e8442b0ab
begin
	histogram(s2, bins = 25, fillalpha = 0.4, title="Estab of diploid migrant")
	xlabel!("\$popsize\$")
	ylabel!("\$n\$")
end

# ╔═╡ 1265d278-030c-48f9-8f2b-b73375b64823
length(s2[s2 .> 0])/1000

# ╔═╡ 84c68f37-fc5a-4383-a142-61aa55df437b
begin
	histogram(s2i, bins = 25, fillalpha = 0.4, title="Estab of diploid migrant")
	xlabel!("\$popsize\$")
	ylabel!("\$n\$")
end

# ╔═╡ 43de9885-feb2-4cac-84ef-1012cef84d87
length(s2i[s2i .> 0])/1000

# ╔═╡ 06a553ee-662d-4bbf-8a5f-365912b14175
island_s4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.)

# ╔═╡ 27568f24-dd03-4f3e-a0cc-296bc9ba2667
island_s4i = IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15., β=0.25)

# ╔═╡ 0d4614ce-3d3b-4b57-ba53-9b42a11c711b
s4 = map(x->(evolving_selectiondeme(MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 0., 0., 1.],1), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.),20)).pop[end], 1:1000) 

# ╔═╡ 288e3007-316f-4c1f-9569-b34e6607e1a4
s4i = map(x->(evolving_islanddeme(IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 0., 0., 1.],1), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15., β=0.25),20)).pop[end], 1:1000) 

# ╔═╡ 345540d4-9ead-4433-922e-24f7c33f7d7f
sort(s4i)

# ╔═╡ 091ec75c-9d54-4b3d-bc45-dbb748ca447a
begin
	histogram(s4, bins = 25, fillalpha = 0.4, title="Estab of tetraploid migrant")
	xlabel!("\$popsize\$")
	ylabel!("\$n\$")
end

# ╔═╡ 558da9c6-e3f1-481e-86c2-8b825ae1f75f
length(s4[s4 .> 0])/1000

# ╔═╡ 637464b6-fe68-4abc-8306-1e784f0352e5
begin
	histogram(s4i, bins = 25, fillalpha = 0.4, title="Estab of tetraploid migrant")
	xlabel!("\$popsize\$")
	ylabel!("\$n\$")
end

# ╔═╡ 45729e96-a64d-46ac-be71-dab33e8e1a9d
length(s4i[s4i .> 0])/1000

# ╔═╡ 71012460-b6af-4282-a78d-692d0a2a9577
md""" IMPORTANT REMARK: The starting phenotypic variance of the island simulated this way for 4N is only half of that of 2N, i.e. there will be less migrants with extreme phenotypes. """

# ╔═╡ 22f70e65-e0b2-4188-91bb-b39c0a74c7a0
md""" ### Single migrant with reduced gamete formation"""

# ╔═╡ 5cdf4ed0-5140-416c-8549-a07cdeaf20e0
island_s2i_UG = IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15., β=0.8)

# ╔═╡ 6abce643-05ba-4635-9f96-c57e981773af
island_p2UG = evolving_islanddeme(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],1), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15., β=0.25), 20)

# ╔═╡ 5c9319b3-158c-4ffc-903a-374ffd4c8cdd
begin
	pf2_p2UG = island_p2UG.p2
	pf3_p2UG = island_p2UG.p3
	pf4_p2UG = island_p2UG.p4
	p21_UG = plot(pf2_p2UG, grid=false, color=:blue, label="diploids", legend=:bottomright)
	#plot!(pf3_p1, grid=false, color=:green, label="triploids")
	plot!(pf4_p2UG, grid=false, color=:red, label="tetraploids")
	hline!([island_s2i_UG.K],label ="K",colour = "black",linestyle=:dash)
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ 7c401941-33bf-4dd5-b231-3ab76c1666d0
begin
	traitmean_ploidy2UG = map(mean, island_p2UG.tm)
	p22_UG = plot(traitmean_ploidy2UG, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright) #title="Stabilizing selection"
	for (i,t) in enumerate(island_p2UG.fta)
	scatter!([i for x in 1:length(t)],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([island_s2i_UG.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ 6da715d6-a372-4962-9d27-e99983783e7c
begin
	p22_fitUG = plot(exp.(island_s2i_UG.β .* (traitmean_ploidy2UG .- island_s2i_UG.θ)), grid=false, color=:red, label=false,linewidth=3,legend=:bottomright) #title="Stabilizing selection"
	for (i,t) in enumerate(island_p2UG.fta)
	scatter!([i for x in 1:length(t)],exp.(island_s2i_UG.β.*(t.-island_s2i_UG.θ)),label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Mean fitness")
	#hline!([island_s2i_UG.θ],label="Optimal henotype",colour="black",linestyle=:dash)
end

# ╔═╡ fd29f45c-f9b8-46ac-afa9-c14456631a32
begin
	Hₒ_ploidyhet2UG = map(mean, 2 .*island_p2UG.het)
	#expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	p23_UG = plot(Hₒ_ploidyhet2UG, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH")
	plot!(1:island_p2UG.ngen, 
		t->expected_heterozygosity(Hₒ_ploidyhet2UG[1], t, 50),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	V_A2UG = map(sum, 2*(0.1)^2 .*island_p2UG.het)
	plot!(V_A2UG, grid=false, color=:red, label="\$V_A(t)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ dbf67254-e59b-4603-9848-e03cdb304f4f
s2i_UG = map(x->(evolving_islanddeme(IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15., β=0.25),20)).pop[end], 1:1000) 

# ╔═╡ 43105e05-0e7e-40f0-b631-6c68ed1b8c3b
begin
	histogram(s2i_UG, bins = 25, fillalpha = 0.4, title="Estab with UG")
	xlabel!("\$popsize\$")
	ylabel!("\$n\$")
end

# ╔═╡ 425c53b6-d7f4-4af6-8051-7bb0baa6cd18
length(s2i_UG[s2i_UG .> 0])/1000

# ╔═╡ 318e6c52-4a13-40af-bc12-af9b1185e6f7
md""" ### Continuous migration """

# ╔═╡ 95ee86cd-14df-42da-b7b6-e094551c298e
function evolving_islandwmigration(d::MixedPloidyDeme, M, ngen; 
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
		migrants = rand(Poisson(M))
		for m in 1:migrants
			migrant = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1)[1]
			push!(d.agents,migrant)
		end
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

# ╔═╡ 4dba1515-8509-402a-97db-437b47cdfda3
function evolving_islanddememigration(d::IslandDeme, M, ngen; 
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
	
	for n=1:ngen
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
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, het=het,tm=tm, af=af, fta=fta) 
end

# ╔═╡ 215b79da-0b74-47f7-8aa2-fb9960165f59
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

# ╔═╡ 5645c34b-5049-49e7-ab47-c0cb33e3e6a8
sim_mig = evolving_islandwmigration(MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5), 1, 20)

# ╔═╡ 6619c52f-5072-42dc-aaa5-1c916e6d040f
sim_migi= evolving_islanddememigration(IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 1, 20)

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

# ╔═╡ 57b440b4-b07b-4875-90e8-7fabae68970e
exp(-1)

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

# ╔═╡ 52129db9-cc6b-4bd9-b3a4-3308a4dfe5f8
begin
	pf2_p2mi = sim_migi.p2
	pf3_p2mi = sim_migi.p3
	pf4_p2mi = sim_migi.p4
	p21mi = plot(pf2_p2mi, grid=false, color=:blue, label="diploids", legend=:topright)
	#plot!(pf3_p1, grid=false, color=:green, label="triploids")
	plot!(pf4_p2mi, grid=false, color=:red, label="tetraploids")
	hline!([island.K],label ="K",colour = "black",linestyle=:dash)
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ bab28060-5865-4693-8ade-ecfd52a3c7c1
begin
	traitmean_migi = map(mean, sim_migi.tm)
	p2mi = plot(traitmean_migi, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright) #title="Stabilizing selection"
	for (i,t) in enumerate(sim_mig.fta)
	scatter!([i for x in 1:length(t)],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([island.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ 72fa64a3-3253-42b9-a8bf-5d47b8fdaa52
begin
	Hₒ_ploidyhet2mi = map(mean, 2 .*sim_migi.het)
	#expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	p23mi = plot(Hₒ_ploidyhet2mi, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH")
	plot!(1:sim_migi.ngen, 
		t->expected_heterozygosity(Hₒ_ploidyhet2mi[1], t, 50),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	V_A2mi = map(sum, 2*(0.1)^2 .*sim_migi.het)
	plot!(V_A2mi, grid=false, color=:red, label="\$V_A(t)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 6efa7a7e-0806-4bde-b519-4753d62fb699
begin
	plot(Poisson(0.01), label="M=0.01", xlim = (0,10),
	ylim = (0,1))
	plot!(Poisson(0.1), label="M=0.1")
	plot!(Poisson(1.), label="M=1")
	plot!(Poisson(10.), label="M=10")
	plot!(Poisson(25.), label="M=25")
	plot!(Poisson(50.), label="M=50")
end

# ╔═╡ 02b516d4-6a80-43d3-8b05-30f56a1459f0
begin
	i_00_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_10_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_20_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_30_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_40_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_50_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_60_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
end

# ╔═╡ 4278d8d1-8634-4be2-be12-0bfe80b67baf
begin
	i_00_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i_10_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i_20_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i_30_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i_40_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25),0.1, 100, 10000)).ngen[end], 1:20)
	i_50_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25),0.1, 100, 10000)).ngen[end], 1:20)
	i_60_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
end

# ╔═╡ e0399c8e-5d34-484d-854b-11fae2c9ee59
begin
	i_00_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_10_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_20_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_30_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_40_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_50_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_60_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
end

# ╔═╡ ce9051fd-2fc7-429a-882d-86fe8f8fc765
begin
	i_00_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i_10_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i_20_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i_30_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i_40_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i_50_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i_60_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
end

# ╔═╡ 3117b134-1168-4318-b845-1bf7f1cbfc4b

# begin
# 	i_00_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i_10_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i_20_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i_30_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i_40_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i_50_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i_60_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# end


# ╔═╡ 96e50f72-39c4-4fa6-9fe3-bb03a202450c
begin
	p2nsim = plot([1,2,3,4,5],[mean(i_00_001b), mean(i_00_01b), mean(i_00_1b), mean(i_00_10b), mean(i_00_100b)], label="△z=0", marker = ([:hex :d]), color=:blue, legend=false)
	plot!([1,2,3,4,5],[mean(i_10_001b), mean(i_10_01b), mean(i_10_1b), mean(i_10_10b), mean(i_10_100b)], label="△z=10", marker = ([:hex :d]), color=:black)
		plot!([1,2,3,4,5],[mean(i_20_001b), mean(i_20_01b), mean(i_20_1b), mean(i_20_10b), mean(i_20_100b)], label="△z=20", marker = ([:hex :d]), color=:blue)
		plot!([1,2,3,4,5],[mean(i_30_001b), mean(i_30_01b), mean(i_30_1b), mean(i_30_10b), mean(i_30_100b)], label="△z=30", marker = ([:hex :d]), color=:black)
		plot!([1,2,3,4,5],[mean(i_40_001b), mean(i_40_01b), mean(i_40_1b), mean(i_40_10b), mean(i_40_100b)], label="△z=40", marker = ([:hex :d]), color=:blue)
		plot!([1,2,3,4,5],[mean(i_50_001b), mean(i_50_01b), mean(i_50_1b), mean(i_50_10b), mean(i_50_100b)], label="△z=50", marker = ([:hex :d]), color=:black)
		plot!([1,2,3,4,5],[mean(i_60_001b), mean(i_60_01b), mean(i_60_1b), mean(i_60_10b), mean(i_60_100b)], label="△z=60", marker = ([:hex :d]), color=:blue)
	xlabel!("Migration rate")
	ylabel!("Time to establish")
#savefig(p2nsim, "Estab_2n_2")
end

# ╔═╡ 04b7b74f-5971-432e-a630-10a67c7f0c3a
begin
	i4_00_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i4_10_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i4_20_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i4_30_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i4_40_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i4_50_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i4_60_1b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
end

# ╔═╡ 9456b032-7a9b-45bb-bfcc-3cad3e23ab7b
begin
	i4_00_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i4_10_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i4_20_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i4_30_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i4_40_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25),0.1, 100, 10000)).ngen[end], 1:20)
	i4_50_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250,[0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25),0.1, 100, 10000)).ngen[end], 1:20)
	i4_60_01b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
end

# ╔═╡ e6f3ec85-47e2-40cb-913f-3c746975b7e4
begin
	i4_00_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i4_10_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i4_20_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i4_30_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i4_40_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i4_50_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i4_60_001b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
end

# ╔═╡ c03cf01a-6741-4d57-a105-3f7bb8a03715
begin
	i4_00_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i4_10_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i4_20_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i4_30_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i4_40_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i4_50_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
	i4_60_10b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 10, 300, 10000)).ngen[end], 1:20)
end

# ╔═╡ 8100b6d0-14d9-4c79-9cfd-d3c9e73fe9bc

# begin
# 	i4_00_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i4_10_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i4_20_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i4_30_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i4_40_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i4_50_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# 	i4_60_100b = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 100, 300, 10000)).ngen[end], 1:20)
# end


# ╔═╡ abf94ddb-4440-4059-b4b4-88b52b0a6f08
begin
	p4nsim = plot([1,2,3,4,5],[mean(i4_00_001b), mean(i4_00_01b), mean(i4_00_1b), mean(i4_00_10b), mean(i4_00_100b)], label="△z=0", marker = ([:hex :d]), color=:blue, legend=false)
	plot!([1,2,3,4,5],[mean(i4_10_001b), mean(i4_10_01b), mean(i4_10_1b), mean(i4_10_10b), mean(i4_10_100b)], label="△z=10", marker = ([:hex :d]), color=:black)
		plot!([1,2,3,4,5],[mean(i4_20_001b), mean(i4_20_01b), mean(i4_20_1b), mean(i4_20_10b), mean(i4_20_100b)], label="△z=20", marker = ([:hex :d]), color=:blue)
		plot!([1,2,3,4,5],[mean(i4_30_001b), mean(i4_30_01b), mean(i4_30_1b), mean(i4_30_10b), mean(i4_30_100b)], label="△z=30", marker = ([:hex :d]), color=:black)
		plot!([1,2,3,4,5],[mean(i4_40_001b), mean(i4_40_01b), mean(i4_40_1b), mean(i4_40_10b), mean(i4_40_100b)], label="△z=40", marker = ([:hex :d]), color=:blue)
		plot!([1,2,3,4,5],[mean(i4_50_001b), mean(i4_50_01b), mean(i4_50_1b), mean(i4_50_10b), mean(i4_50_100b)], label="△z=50", marker = ([:hex :d]), color=:black)
	plot!([1,2,3,4,5],[mean(i4_60_001b), mean(i4_60_01b), mean(i4_60_1b), mean(i4_60_10b), mean(i4_60_100b)], label="△z=50", marker = ([:hex :d]), color=:blue)
	xlabel!("Migration rate")
	ylabel!("Time to establish")
#savefig(p4nsim, "Estab_4n_2")
end

# ╔═╡ 12bd4301-7a63-482a-be63-bf5f7c78bb90
md""" ### Cytotype load"""

# ╔═╡ 5452b1f1-78f9-4ed4-9050-c2867881872f
u = 0.1

# ╔═╡ 599b73d6-b3df-4e0e-b227-1594f062ecc1
begin
	i_00_1bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_10_1bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_20_1bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_30_1bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_40_1bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_50_1bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
	i_60_1bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 1, 100, 10000)).ngen[end], 1:20)
end

# ╔═╡ da88070d-da7a-45b7-a756-f85bd1a63dda
begin
	i_00_01bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i_10_01bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i_20_01bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i_30_01bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i_40_01bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i_50_01bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
	i_60_01bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 0.1, 100, 10000)).ngen[end], 1:20)
end

# ╔═╡ ce586b06-6a17-4360-8ab1-be60b90087a0
begin
	i_00_001bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_10_001bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_20_001bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_30_001bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_40_001bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_50_001bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
	i_60_001bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 0.01, 100, 10000)).ngen[end], 1:20)
end

# ╔═╡ 70f024ed-0a5a-4630-a087-1e2ae0c9627e
begin
	i_00_10bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 10., 100, 10000)).ngen[end], 1:20)
	i_10_10bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 10., 100, 10000)).ngen[end], 1:20)
	i_20_10bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 10., 100, 10000)).ngen[end], 1:20)
	i_30_10bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 10, 100, 10000)).ngen[end], 1:20)
	i_40_10bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 10., 100, 10000)).ngen[end], 1:20)
	i_50_10bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 10., 100, 10000)).ngen[end], 1:20)
	i_60_10bUG = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 10., 100, 10000)).ngen[end], 1:20)
end

# ╔═╡ 15a60434-8fac-4f55-a98a-958de300ce3c
begin
	p2nsimUG = plot([1,2,3,4],[mean(i_00_001bUG), mean(i_00_01bUG), mean(i_00_1bUG), mean(i_00_10bUG)], label="△z=0", marker = ([:hex :d]), color=:blue, title="u=0.1")
	plot!([1,2,3,4],[mean(i_10_001bUG), mean(i_10_01bUG), mean(i_10_1bUG), mean(i_10_10bUG)], label="△z=1", marker = ([:hex :d]), color=:black)
		plot!([1,2,3,4],[mean(i_20_001bUG), mean(i_20_01bUG), mean(i_20_1bUG), mean(i_20_10bUG)], label="△z=2", marker = ([:hex :d]), color=:blue)
		plot!([1,2,3,4],[mean(i_30_001bUG), mean(i_30_01bUG), mean(i_30_1bUG), mean(i_30_10bUG)], label="△z=3", marker = ([:hex :d]), color=:black)
		plot!([1,2,3,4],[mean(i_40_001bUG), mean(i_40_01bUG), mean(i_40_1bUG), mean(i_40_10bUG)], label="△z=4", marker = ([:hex :d]), color=:blue)
		plot!([1,2,3,4],[mean(i_50_001bUG), mean(i_50_01bUG), mean(i_50_1bUG), mean(i_50_10bUG)], label="△z=5", marker = ([:hex :d]), color=:black)
	plot!([1,2,3,4],[mean(i_60_001bUG), mean(i_60_01bUG), mean(i_60_1bUG), mean(i_60_10bUG)], label="△z=6", marker = ([:hex :d]), color=:blue)
	xlabel!("Migration rate")
	ylabel!("Time to establish")
#savefig(p2nsimUG, "Estab_2n_UG020020")
end

# ╔═╡ 093d6970-741f-4a02-a15d-e7fa6d44833a
begin
	i_00_1bUG2 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 1., 100, 10000)).p2[end], 1:20)
	i_10_1bUG2 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 1., 100, 10000)).p2[end], 1:20)
	i_20_1bUG2 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 1., 100, 10000)).p2[end], 1:20)
	i_30_1bUG2 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 1., 100, 10000)).p2[end], 1:20)
	i_40_1bUG2 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 1., 100, 10000)).p2[end], 1:20)
	i_50_1bUG2 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 1., 100, 10000)).p2[end], 1:20)
	i_60_1bUG2 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 1., 100, 10000)).p2[end], 1:20)
end

# ╔═╡ 94bc5982-14b2-406a-aa4d-9fcc41cbd7c2
begin
	i_00_1bUG4 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=12.5, β=0.25), 1., 100, 10000)).p4[end], 1:20)
	i_10_1bUG4 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=13.5, β=0.25), 1., 100, 10000)).p4[end], 1:20)
	i_20_1bUG4 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=14.5, β=0.25), 1., 100, 10000)).p4[end], 1:20)
	i_30_1bUG4 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15.5, β=0.25), 1., 100, 10000)).p4[end], 1:20)
	i_40_1bUG4 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=16.5, β=0.25), 1., 100, 10000)).p4[end], 1:20)
	i_50_1bUG4 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=17.5, β=0.25), 1., 100, 10000)).p4[end], 1:20)
	i_60_1bUG4 = map(x->(evolving_islandwbreak(IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=18.5, β=0.25), 1., 100, 10000)).p4[end], 1:20)
end

# ╔═╡ 55776bc3-55a3-4890-bb4c-49e988a26da3
begin
	p2nsimUGp24 = plot([1,2,3,4,5,6,7],[mean(i_00_1bUG2), mean(i_10_1bUG2), mean(i_20_1bUG2), mean(i_30_1bUG2), mean(i_40_1bUG2), mean(i_50_1bUG2), mean(i_60_1bUG2)], marker = ([:hex :d]), color=:blue, label="diploid", title="u=0.1, m=1")
	
	plot!([1,2,3,4,5,6,7],[mean(i_00_1bUG4), mean(i_10_1bUG4), mean(i_20_1bUG4), mean(i_30_1bUG4), mean(i_40_1bUG4), mean(i_50_1bUG4), mean(i_60_1bUG4)], marker = ([:hex :d]), color=:black, label="tetraploid")

	xlabel!("△z")
	ylabel!("Population size")
#savefig(p2nsimUG, "Estab_2n_UG020020")
end

# ╔═╡ ae97a439-2eb8-43c8-a983-8cc4eaff3d9b
function grid_search(t)
	ploidy = []
	param = []
	pop_size = []
	p2 = []
	p4 = []
	for u in range(0, stop=0.5, length=t)
		for rep in 1:20
		UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.]
		d_p = IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UG, μ=0., θ=14.5, β=0.25)
		sim_ploidyvar = evolving_islandwbreak(d_p, .01, 100, 10000)
		if sim_ploidyvar.p2[end] >= sim_ploidyvar.p4[end]
			push!(ploidy,2)
		else
			push!(ploidy,4)
		end
		push!(param, u)
		pop = (sim_ploidyvar.p2[end] + sim_ploidyvar.p4[end])
		push!(pop_size, pop)
		push!(p2, sim_ploidyvar.p2[end])
		push!(p4, sim_ploidyvar.p4[end])
		end
	end
	ploidy, param, pop_size, p2, p4
end

# ╔═╡ b330afb1-a2be-4f6d-9395-93b2ff931031
stats = grid_search(200)

# ╔═╡ ab1bb2e8-2317-45b0-a62d-b640526d96ff
dp_2 = [(stats[2][i],stats[1][i],stats[3][i]) for i in 1:2000]

# ╔═╡ 13f6f773-5dd4-421b-9b35-7ae09bd0c3d7
begin
function prob(b)
		c = 0
		for x in b
			if x == 4
				c += 1
			end
		end
		c/length(b)
	end

function stabprob(a)
	i = 1
	j = 20
	p = []
	while j <= length(a)
		push!(p,prob(a[i:j]))
		i += 20
		j += 20
	end
	p
	end	
end	

# ╔═╡ 810211aa-0fbf-41e9-b0c8-c01fce34adcb
begin
tick1 = stabprob(stats[1])
plot([0.0025:0.0025:0.5...],tick1,label=false)
vline!([0.17],label="u=0.17",linewidth=2,style=:dash)
hline!([0.50],label=false,linewidth=2,style=:dash)
xlabel!("u")
ylabel!("P estab")
end

# ╔═╡ 43560e72-ad13-4efe-9f03-30d19a2ddf25
begin
p3ug = scatter(stats[2],stats[1], label=false, title="UG for 2n and 4n")
vline!([0.17], label="u=0.17",linewidth=5)
xlabel!("u (unreduced gamete formation 2n)")
ylabel!("Ploidy")
end

# ╔═╡ 3a8cbad9-f3d6-4e34-afd2-783774c00af1
function bin(data)
	up2 = []
	up4 = []
	for i in data
		if i[3] != 0
			if i[2] == 2
				push!(up2,i[1])
			else
				push!(up4,i[1])
			end
		end
	end
	up2,up4
end

# ╔═╡ 69b3670a-e159-4434-88c1-4142a36eec6f
binz = bin(dp_2)

# ╔═╡ c4382312-f4ce-40de-857a-4c25122c6098
begin
	p1ug = histogram(binz[1], bins=50, label=false, title="Diploid > Tetraploid")
	vline!([0.17], label="u=0.17",linewidth=5)
	xlabel!("u")
	ylabel!("Count")
end

# ╔═╡ e552c371-4288-433f-a3f3-993cc6bd8e35
begin
	p2ug = histogram(binz[2], bins=50, label=false, title="Tetraploid > Diploid")
	vline!([0.17], label="u=0.17",linewidth=5)
	xlabel!("u")
	ylabel!("Count")
end

# ╔═╡ 2d4788a6-b93b-4b3d-87b5-354d723f7499
begin
p7 = scatter(stats[2], stats[3], grid=false, color=:black, label="Pop size after t generations")
scatter!(stats[2], stats[4], grid=false, color=:green, label="Diploids")
scatter!(stats[2], stats[5], grid=false, color=:red, label="Tetraploids")
vline!([0.17], label="u=0.17",linewidth=5)
xlabel!("\$u\$")
ylabel!("Number of individuals")
end

# ╔═╡ 8358bab3-2640-4830-88e0-ba0413333551
md""" Is this only for pops that estab >100? effect + effect of dir sel vs stab sel"""

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
# ╠═bfc7fa55-15a7-4690-9c33-7cfffedc975c
# ╠═6c24a8d9-7ecd-471a-b3d8-7b7ab6ef32d2
# ╠═09247474-d2d7-4a19-b3b6-dccf77fcd7f0
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
# ╟─4ce3fdbb-82f2-46da-a226-485306b0d91b
# ╟─e87406ec-2deb-4fda-9e49-5f95676b32c3
# ╟─602cebcc-46bf-48fe-a3d3-163318796b98
# ╠═1182e78a-f422-4e14-b6f4-e6f2c05c3da4
# ╠═824018c1-6f6b-40ff-8b77-aebaccc20afc
# ╠═8841a695-1110-487e-bd8a-720ef80a64c3
# ╠═2da9fafc-64c4-4c4c-a0c8-2e248ae84441
# ╠═ef4e8ee0-828f-478c-ba12-1a555d620019
# ╠═b3f9e34a-2a82-4223-ac2a-6c9e8442b0ab
# ╠═1265d278-030c-48f9-8f2b-b73375b64823
# ╠═84c68f37-fc5a-4383-a142-61aa55df437b
# ╠═43de9885-feb2-4cac-84ef-1012cef84d87
# ╠═06a553ee-662d-4bbf-8a5f-365912b14175
# ╠═27568f24-dd03-4f3e-a0cc-296bc9ba2667
# ╠═0d4614ce-3d3b-4b57-ba53-9b42a11c711b
# ╠═288e3007-316f-4c1f-9569-b34e6607e1a4
# ╠═345540d4-9ead-4433-922e-24f7c33f7d7f
# ╠═091ec75c-9d54-4b3d-bc45-dbb748ca447a
# ╠═558da9c6-e3f1-481e-86c2-8b825ae1f75f
# ╠═637464b6-fe68-4abc-8306-1e784f0352e5
# ╠═45729e96-a64d-46ac-be71-dab33e8e1a9d
# ╟─71012460-b6af-4282-a78d-692d0a2a9577
# ╟─22f70e65-e0b2-4188-91bb-b39c0a74c7a0
# ╠═5cdf4ed0-5140-416c-8549-a07cdeaf20e0
# ╠═6abce643-05ba-4635-9f96-c57e981773af
# ╠═5c9319b3-158c-4ffc-903a-374ffd4c8cdd
# ╠═7c401941-33bf-4dd5-b231-3ab76c1666d0
# ╠═6da715d6-a372-4962-9d27-e99983783e7c
# ╠═fd29f45c-f9b8-46ac-afa9-c14456631a32
# ╠═dbf67254-e59b-4603-9848-e03cdb304f4f
# ╠═43105e05-0e7e-40f0-b631-6c68ed1b8c3b
# ╠═425c53b6-d7f4-4af6-8051-7bb0baa6cd18
# ╟─318e6c52-4a13-40af-bc12-af9b1185e6f7
# ╠═44771223-a589-4de3-9fea-321b70138d80
# ╠═c0a08e7f-dca4-4107-b2b3-d18d7dae77a8
# ╠═95ee86cd-14df-42da-b7b6-e094551c298e
# ╠═4dba1515-8509-402a-97db-437b47cdfda3
# ╠═215b79da-0b74-47f7-8aa2-fb9960165f59
# ╠═5645c34b-5049-49e7-ab47-c0cb33e3e6a8
# ╠═6619c52f-5072-42dc-aaa5-1c916e6d040f
# ╠═f1595c57-e039-41aa-b8e9-9997cdc19641
# ╠═57b440b4-b07b-4875-90e8-7fabae68970e
# ╠═ce63121d-7c50-49d2-a60c-1771bcd01ae9
# ╠═f51ad793-1dc2-4515-8b7c-627904bc016e
# ╟─52129db9-cc6b-4bd9-b3a4-3308a4dfe5f8
# ╟─bab28060-5865-4693-8ade-ecfd52a3c7c1
# ╠═72fa64a3-3253-42b9-a8bf-5d47b8fdaa52
# ╠═6efa7a7e-0806-4bde-b519-4753d62fb699
# ╟─02b516d4-6a80-43d3-8b05-30f56a1459f0
# ╟─4278d8d1-8634-4be2-be12-0bfe80b67baf
# ╟─e0399c8e-5d34-484d-854b-11fae2c9ee59
# ╟─ce9051fd-2fc7-429a-882d-86fe8f8fc765
# ╟─3117b134-1168-4318-b845-1bf7f1cbfc4b
# ╟─96e50f72-39c4-4fa6-9fe3-bb03a202450c
# ╟─04b7b74f-5971-432e-a630-10a67c7f0c3a
# ╟─9456b032-7a9b-45bb-bfcc-3cad3e23ab7b
# ╟─e6f3ec85-47e2-40cb-913f-3c746975b7e4
# ╟─c03cf01a-6741-4d57-a105-3f7bb8a03715
# ╟─8100b6d0-14d9-4c79-9cfd-d3c9e73fe9bc
# ╟─abf94ddb-4440-4059-b4b4-88b52b0a6f08
# ╟─12bd4301-7a63-482a-be63-bf5f7c78bb90
# ╠═5452b1f1-78f9-4ed4-9050-c2867881872f
# ╠═599b73d6-b3df-4e0e-b227-1594f062ecc1
# ╠═da88070d-da7a-45b7-a756-f85bd1a63dda
# ╠═ce586b06-6a17-4360-8ab1-be60b90087a0
# ╠═70f024ed-0a5a-4630-a087-1e2ae0c9627e
# ╠═15a60434-8fac-4f55-a98a-958de300ce3c
# ╠═093d6970-741f-4a02-a15d-e7fa6d44833a
# ╠═94bc5982-14b2-406a-aa4d-9fcc41cbd7c2
# ╠═55776bc3-55a3-4890-bb4c-49e988a26da3
# ╠═ae97a439-2eb8-43c8-a983-8cc4eaff3d9b
# ╠═b330afb1-a2be-4f6d-9395-93b2ff931031
# ╠═ab1bb2e8-2317-45b0-a62d-b640526d96ff
# ╠═13f6f773-5dd4-421b-9b35-7ae09bd0c3d7
# ╠═810211aa-0fbf-41e9-b0c8-c01fce34adcb
# ╠═43560e72-ad13-4efe-9f03-30d19a2ddf25
# ╠═3a8cbad9-f3d6-4e34-afd2-783774c00af1
# ╠═69b3670a-e159-4434-88c1-4142a36eec6f
# ╠═c4382312-f4ce-40de-857a-4c25122c6098
# ╠═e552c371-4288-433f-a3f3-993cc6bd8e35
# ╠═2d4788a6-b93b-4b3d-87b5-354d723f7499
# ╟─8358bab3-2640-4830-88e0-ba0413333551
