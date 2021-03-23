### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 1f791df0-8b5e-11eb-3bc5-e7302263af09
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes

# ╔═╡ 325d63e0-8b5e-11eb-0d55-f149adc308e1
using PolyStab

# ╔═╡ 423cffe0-8b5f-11eb-1df4-ff7eed10289a
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_deme_UG, evolving_deme_ploidyvar, heterozygosities_p, allelefreqs_p, random_mating_mixedp, neutral_evolving_deme

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
island = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 0.1, 250, [2], 0, d=2.),randagent_p(0.5, 0.1, 250, [4], 0, d=4.)), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], K=100 )

# ╔═╡ 8272fcde-8b5f-11eb-1c13-b591229f902f
migrant = randagent_p(0.48, 0.1, 250, [2], 1, d=2.)[1]

# ╔═╡ de4112e2-8b60-11eb-1e75-1ffaf3a0c06d
trait(migrant)

# ╔═╡ faa45f90-8b6b-11eb-2826-776161a15a58
migrant.loci

# ╔═╡ ce694600-8bd6-11eb-0d50-c7bd17294608
allelefreqs_p(migrant)

# ╔═╡ 05a7f6f0-8bd9-11eb-0892-f3cde611f9ee
mean(allelefreqs_p(migrant))

# ╔═╡ 143183d0-8b61-11eb-1986-310271688189
push!(island.agents,migrant)

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

# ╔═╡ 58605950-8b61-11eb-16d9-0ff9be46d489
island_p = evolving_deme_ploidyvar(island,100)

# ╔═╡ bc1c99b0-8b69-11eb-29f7-83efe7ac0372
island_p.deme

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

# ╔═╡ ee422a80-8b6a-11eb-13ce-a7f8ab3de0ef
begin
	Hₒ_ploidyhet = map(mean, 2 .*island_p.het)
	expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
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

# ╔═╡ 84585d6e-8c09-11eb-3c5c-af38eb556b15
var_add(a::Agent,α) = ploidy(a)*α^2*sum(heterozygosities_p(a))

# ╔═╡ 5a0f4c50-8c26-11eb-0f10-dfd79a2ff7c5
begin
	plot(p1,p2,p3)
	savefig("single_migrant.pdf")
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
# ╠═8272fcde-8b5f-11eb-1c13-b591229f902f
# ╠═de4112e2-8b60-11eb-1e75-1ffaf3a0c06d
# ╠═faa45f90-8b6b-11eb-2826-776161a15a58
# ╠═ce694600-8bd6-11eb-0d50-c7bd17294608
# ╠═05a7f6f0-8bd9-11eb-0892-f3cde611f9ee
# ╠═143183d0-8b61-11eb-1986-310271688189
# ╠═4d1ba2c0-8bcf-11eb-2b00-adaff59d3e9e
# ╠═7b190320-8bcf-11eb-15b6-2b3a871bc3bc
# ╠═93164980-8bcd-11eb-02a1-496109bc20c3
# ╠═91834c70-8bd8-11eb-0ef6-e79beaa1849b
# ╠═52c76560-8b61-11eb-1edd-e71d92ca9a1e
# ╠═58605950-8b61-11eb-16d9-0ff9be46d489
# ╠═bc1c99b0-8b69-11eb-29f7-83efe7ac0372
# ╠═6df662b0-8b6a-11eb-0e74-611b7747413e
# ╠═dc3ad710-8b6a-11eb-3e55-7ba4760f8132
# ╠═ee422a80-8b6a-11eb-13ce-a7f8ab3de0ef
# ╠═84585d6e-8c09-11eb-3c5c-af38eb556b15
# ╠═5a0f4c50-8c26-11eb-0f10-dfd79a2ff7c5
