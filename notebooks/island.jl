### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 1f791df0-8b5e-11eb-3bc5-e7302263af09
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes

# ╔═╡ 325d63e0-8b5e-11eb-0d55-f149adc308e1
using PolyStab

# ╔═╡ 423cffe0-8b5f-11eb-1df4-ff7eed10289a
using PolyStab:randagent_p, MixedPloidyDeme, trait, evolving_deme_UG, evolving_deme_ploidyvar

# ╔═╡ 7a5c3050-8b62-11eb-1d96-3da8035911e4
using StatsBase:weights

# ╔═╡ a7bc4e02-8b5c-11eb-0e36-ef47de093003
md""" ### Island model
"""

# ╔═╡ ba38834e-8b5c-11eb-246c-99429ab2ae27
md""" ###### cfr. Establishment in a new habitat by polygenic adaptation (Barton and Etheridge, 2017)"""

# ╔═╡ baffcaf2-8b5c-11eb-1479-e110e25be9d5
md""" Representation: large abstract mainland population (assumes: infinite population size in HWE), migration from mainland to different smaller islands with finite population size (respresented as demes).

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
island = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 0.1, 250, [2], 0, d=2.),randagent_p(0.5, 0.1, 250, [4], 0, d=4.)), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0.9 0.1 0. 0. ; 0. 0. 0. 0. ; 0. 0.9 0. 0.1] )

# ╔═╡ 8272fcde-8b5f-11eb-1c13-b591229f902f
migrant = randagent_p(0.45, 0.1, 250, [2], 1, d=2.)[1]

# ╔═╡ faa45f90-8b6b-11eb-2826-776161a15a58


# ╔═╡ de4112e2-8b60-11eb-1e75-1ffaf3a0c06d
trait(migrant)

# ╔═╡ 143183d0-8b61-11eb-1986-310271688189
push!(island.agents,migrant)

# ╔═╡ 52c76560-8b61-11eb-1edd-e71d92ca9a1e
island

# ╔═╡ 58605950-8b61-11eb-16d9-0ff9be46d489
island_p = evolving_deme_ploidyvar(island,50)

# ╔═╡ bc1c99b0-8b69-11eb-29f7-83efe7ac0372
island_p.deme

# ╔═╡ 6df662b0-8b6a-11eb-0e74-611b7747413e
begin
	pf2_p1 = island_p.p2
	pf3_p1 = island_p.p3
	pf4_p1 = island_p.p4
	plot(pf2_p1, grid=false, color=:blue, label="diploid")
	plot!(pf3_p1, grid=false, color=:green, label="triploid")
	plot!(pf4_p1, grid=false, color=:red, label="tetraploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ dc3ad710-8b6a-11eb-3e55-7ba4760f8132
begin
	traitmean_ploidy = map(mean, island_p.tm)
	plot(traitmean_ploidy, grid=false, color=:black, label=false, title="Effect of stabilizing selection")
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ ee422a80-8b6a-11eb-13ce-a7f8ab3de0ef
begin
	Hₒ_ploidyvar = map(mean, island_p.het)
	expected_heterozygosity(H₀, t, N) = ((1.0-1.0/N)^t)*H₀
	plot(Hₒ_ploidyvar, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH mixed ploidy with selection")
	plot!(0:island_p.ngen+1, 
		t->expected_heterozygosity(2*0.5*(0.5), t, 100),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/2N)^t H_o(0)\$")
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
# ╠═7a5c3050-8b62-11eb-1d96-3da8035911e4
# ╟─3712a3fe-8b5d-11eb-2560-5940bd21ea0c
# ╟─d97a6a70-8b5d-11eb-0ed6-056d834ce9bb
# ╠═2625fd72-8b5f-11eb-1e69-adf12cbd84c9
# ╠═8272fcde-8b5f-11eb-1c13-b591229f902f
# ╠═faa45f90-8b6b-11eb-2826-776161a15a58
# ╠═de4112e2-8b60-11eb-1e75-1ffaf3a0c06d
# ╠═143183d0-8b61-11eb-1986-310271688189
# ╠═52c76560-8b61-11eb-1edd-e71d92ca9a1e
# ╠═58605950-8b61-11eb-16d9-0ff9be46d489
# ╠═bc1c99b0-8b69-11eb-29f7-83efe7ac0372
# ╠═6df662b0-8b6a-11eb-0e74-611b7747413e
# ╠═dc3ad710-8b6a-11eb-3e55-7ba4760f8132
# ╠═ee422a80-8b6a-11eb-13ce-a7f8ab3de0ef
