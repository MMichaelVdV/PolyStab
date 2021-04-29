### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ b1ccc66f-880d-40c5-8309-84e31c2fe1d1
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, DataFrames, StatsPlots

# ╔═╡ 81e0696e-1d12-410c-89e8-86e492810a62
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, evolving_neutraldeme, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy, trait_mean, randagent, evolving_haploiddeme, mate_p, evolving_selectiondemeh

# ╔═╡ 588b9cb0-a929-11eb-2ac2-e77ecad3b43f
md""" ##### Some notes on the Fundamental Theorem of Natural Selection. """

# ╔═╡ d2107ebd-72c0-48fb-8dbb-7a6478780104
md""" Phenotypic variance seems to persist longer in tetraploids under stabilizing selection, this coulbd be mainly due to the effect of less genetic drift or something else. Tetraploid population takes longer to reach equilibrium."""

# ╔═╡ 1a688141-c921-44c1-9685-818aff3c7b9e
md"""The rate of return to equilibrium population size in PnB (2015) is equivalent to r_m - V_g / 2 V_s ."""

# ╔═╡ c1d263e3-d3d7-48f4-871f-05f82f858b8b
begin
d_p1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [1., 0., 0., 0.],150), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 0.1)
d_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [0., 1., 0., 0.],150), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 0.1)
d_p4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [0., 0., 0., 1.],150), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 0.1)
end

# ╔═╡ c1c22e59-3fb6-45df-b793-3aa463398aec
begin
stabsel_p1 = evolving_selectiondemeh(d_p1,500)
stabsel_p2 = evolving_selectiondeme(d_p2,500)
stabsel_p4 = evolving_selectiondeme(d_p4,500)
end

# ╔═╡ 9e36a5a1-6461-4e8a-b088-cc838f86dd5e
begin
plot(Normal(mean(stabsel_p1.fta[1]), std(stabsel_p1.fta[1])), color=:black, label="Haploid")
plot!(Normal(mean(stabsel_p2.fta[1]), std(stabsel_p2.fta[1])), color=:blue, label="Diploid")
plot!(Normal(mean(stabsel_p4.fta[1]), std(stabsel_p4.fta[1])), color=:red, label="Tetraploid")
xlabel!("Phenotype before selection")
end

# ╔═╡ 898ea668-2254-4927-8fcb-fdb1ccf3e6ec
begin
	popsize_p1 = map(mean, stabsel_p1.pop)
	popsize_p2 = map(mean, stabsel_p2.pop)
	popsize_p4 = map(mean, stabsel_p4.pop)
	plot(popsize_p1, grid=false, color=:black, label=false)
	plot!(popsize_p2, grid=false, color=:blue, label=false)
	plot!(popsize_p4, grid=false, color=:red, label=false)
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ a5f6483f-1567-4ebf-b58e-a08b229b2ef9
begin
	popp1 = plot(popsize_p1, grid=false, color=:black, label=false)
	popp2 = plot(popsize_p2, grid=false, color=:blue, label=false)
	popp4 = plot(popsize_p4, grid=false, color=:red, label=false)
	plot(popp1,popp2,popp4, layout = (3, 1))
end

# ╔═╡ 1e6a7e16-3585-4135-a0fc-20e9c89105ea
md""" #### Functions """

# ╔═╡ Cell order:
# ╟─588b9cb0-a929-11eb-2ac2-e77ecad3b43f
# ╟─d2107ebd-72c0-48fb-8dbb-7a6478780104
# ╠═1a688141-c921-44c1-9685-818aff3c7b9e
# ╠═c1d263e3-d3d7-48f4-871f-05f82f858b8b
# ╠═c1c22e59-3fb6-45df-b793-3aa463398aec
# ╠═9e36a5a1-6461-4e8a-b088-cc838f86dd5e
# ╠═898ea668-2254-4927-8fcb-fdb1ccf3e6ec
# ╠═a5f6483f-1567-4ebf-b58e-a08b229b2ef9
# ╟─1e6a7e16-3585-4135-a0fc-20e9c89105ea
# ╠═b1ccc66f-880d-40c5-8309-84e31c2fe1d1
# ╠═81e0696e-1d12-410c-89e8-86e492810a62
