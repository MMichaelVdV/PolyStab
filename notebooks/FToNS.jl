### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ b1ccc66f-880d-40c5-8309-84e31c2fe1d1
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, DataFrames, StatsPlots

# ╔═╡ 81e0696e-1d12-410c-89e8-86e492810a62
using PolyStab: Agent, randagent_p, MixedPloidyDeme, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, evolving_neutraldeme, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy, trait_mean, randagent, evolving_haploiddeme, mate_p, evolving_selectiondemeh

# ╔═╡ 588b9cb0-a929-11eb-2ac2-e77ecad3b43f
md""" ##### Some notes on genetic load and stabilizing selection. """

# ╔═╡ d2107ebd-72c0-48fb-8dbb-7a6478780104
md""" Phenotypic variance seems to persist longer in tetraploids under stabilizing selection, this coulbd be mainly due to the effect of less genetic drift or something else. Tetraploid population takes longer to reach equilibrium."""

# ╔═╡ 1a688141-c921-44c1-9685-818aff3c7b9e
md"""The rate of return to equilibrium population size in PnB (2015) is equivalent to r_m - V_g / 2 V_s ."""

# ╔═╡ c1d263e3-d3d7-48f4-871f-05f82f858b8b
begin
d_p1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [1., 0., 0., 0.],150), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)
d_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [0., 1., 0., 0.],150), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)
d_p4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [0., 0., 0., 1.],150), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)
end

# ╔═╡ 8fa49926-3bf2-4eb9-95c8-634941600ea2
begin
	
d25_p1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [1., 0., 0., 0.],150, d=0.25), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)
d25_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [0., 1., 0., 0.],150, d=0.25), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)
d25_p4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [0., 0., 0., 1.],150, d=0.25), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)

d50_p1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [1., 0., 0., 0.],150, d=0.50), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)
d50_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [0., 1., 0., 0.],150, d=0.50), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)
d50_p4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [0., 0., 0., 1.],150, d=0.50), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)

d75_p1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [1., 0., 0., 0.],150, d=0.75), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)
d75_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [0., 1., 0., 0.],150, d=0.75), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)
d75_p4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [0., 0., 0., 1.],150, d=0.75), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = .25, α=0.8)

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

# ╔═╡ e7b41fc7-af1a-4b0e-aff7-f21bb85d1e9f
md""" mean pop size of tetraploids stays lower? """

# ╔═╡ a5f6483f-1567-4ebf-b58e-a08b229b2ef9
begin
	popp1 = plot(popsize_p1, grid=false, color=:black, label=false)
	popp2 = plot(popsize_p2, grid=false, color=:blue, label=false)
	popp4 = plot(popsize_p4, grid=false, color=:red, label=false)
	plot(popp1,popp2,popp4, layout = (3, 1))
end

# ╔═╡ e5596486-eab1-4a86-a31d-c33ded507add
md""" #### Effect of allelic effect size 'α' and number of loci 'n' """

# ╔═╡ 3098bc8d-42a5-4ff9-9931-9d7553a3b489
md""" but also effect of other parameters: Vs, rm. """

# ╔═╡ ae903443-4310-45ce-beaa-53aec00db7ab
begin
d_p2_1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.1, 400, [0., 1., 0., 0.],250), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 2., K=500, α=0.1)
d_p2_2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [0., 1., 0., 0.],250), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 2., K=500, α=0.2)
d_p2_3 = MixedPloidyDeme(agents = randagent_p(0.5, 0.4, 100, [0., 1., 0., 0.],250), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 2., K=500, α=0.4)
d_p2_4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.8, 50, [0., 1., 0., 0.],250), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 2., K=500, α=0.8)
end

# ╔═╡ 5f227e3d-581c-4b50-aad9-92931427e3b2
begin
stabsel_p2_1 = evolving_selectiondeme(d_p2_1,500)
stabsel_p2_2 = evolving_selectiondeme(d_p2_2,500)
stabsel_p2_3 = evolving_selectiondeme(d_p2_3,500)
stabsel_p2_4 = evolving_selectiondeme(d_p2_4,500)
end

# ╔═╡ 32f5be62-7fb2-400f-961e-8b394d91cb97
begin
plot(Normal(mean(stabsel_p2_1.fta[1]), std(stabsel_p2_1.fta[1])), color=:black, label="α=0.1, n=400")
plot!(Normal(mean(stabsel_p2_2.fta[1]), std(stabsel_p2_2.fta[1])), color=:blue, label="α=0.2, n=200")
plot!(Normal(mean(stabsel_p2_3.fta[1]), std(stabsel_p2_3.fta[1])), color=:red, label="α=0.4, n=100")
plot!(Normal(mean(stabsel_p2_4.fta[1]), std(stabsel_p2_4.fta[1])), color=:green, label="α=0.8, n=50")
xlabel!("Phenotype before selection")
end

# ╔═╡ e549e9e6-2d83-4fdb-967b-b668cc209687
begin
	popsize_p2_1 = map(mean, stabsel_p2_1.pop)
	popsize_p2_2 = map(mean, stabsel_p2_2.pop)
	popsize_p2_3 = map(mean, stabsel_p2_3.pop)
	popsize_p2_4 = map(mean, stabsel_p2_4.pop)
	plot(popsize_p2_1, grid=false, color=:black, label=false)
	plot!(popsize_p2_2, grid=false, color=:blue, label=false)
	plot!(popsize_p2_3, grid=false, color=:red, label=false)
	plot!(popsize_p2_4, grid=false, color=:green, label=false)
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ 5b22b365-1898-4079-bd53-a164ed2f94fc
begin
	popp2_1 = plot(popsize_p2_1, grid=false, color=:black, label=false)
	popp2_2 = plot(popsize_p2_2, grid=false, color=:blue, label=false)
	popp2_3 = plot(popsize_p2_3, grid=false, color=:red, label=false)
	popp2_4 = plot(popsize_p2_4, grid=false, color=:green, label=false)
	plot(popp2_1,popp2_2,popp2_3,popp2_4, layout = (4, 1))
end

# ╔═╡ 9b6e587e-6266-42b8-853e-137d135422bd
md""" Why do these populations crash to 0 ? Some effect of stabilizing selection and initialization."""

# ╔═╡ ded2e881-c5a2-4faa-83c0-0471e633a522
begin
	Hₒs_p1 = map(mean, stabsel_p2_1.het)
	Hₒs_p2 = map(mean, stabsel_p2_2.het)
	Hₒs_p3 = map(mean, stabsel_p2_3.het)
	Hₒs_p4 = map(mean, stabsel_p2_4.het)

	plot(Hₒs_p1, grid=false, color=:black, label="α=0.1, n=400", title="LOH")
	plot!(Hₒs_p2, grid=false, color=:blue, label="α=0.2, n=200")
	plot!(Hₒs_p3, grid=false, color=:red, label="α=0.4, n=100")
	plot!(Hₒs_p4, grid=false, color=:green, label="α=0.8, n=50")

	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 7a021230-ccab-43bd-bef1-5ce25e9db559
begin
neutral_p2_1 = evolving_neutraldeme(d_p2_1,500)
neutral_p2_2 = evolving_neutraldeme(d_p2_2,500)
neutral_p2_3 = evolving_neutraldeme(d_p2_3,500)
neutral_p2_4 = evolving_neutraldeme(d_p2_4,500)
end

# ╔═╡ 5b413f3d-f2dd-48fc-b761-1528b4604d7a
begin
	Hₒn_p1 = map(mean, neutral_p2_1.het)
	Hₒn_p2 = map(mean, neutral_p2_2.het)
	Hₒn_p3 = map(mean, neutral_p2_3.het)
	Hₒn_p4 = map(mean, neutral_p2_4.het)

	plot(Hₒn_p1, grid=false, color=:black, label="α=0.1, n=400", title="LOH")
	plot!(Hₒn_p2, grid=false, color=:blue, label="α=0.2, n=200")
	plot!(Hₒn_p3, grid=false, color=:red, label="α=0.4, n=100")
	plot!(Hₒn_p4, grid=false, color=:green, label="α=0.8, n=50")

	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 1e6a7e16-3585-4135-a0fc-20e9c89105ea
md""" #### Functions """

# ╔═╡ 4bd0cf27-361a-4ec8-811e-f20f019c9c15
"""
    trait(a::Agent)
"""
trait(a::Agent) = a.d*sum(a)/ploidy(a)

# ╔═╡ 367df54e-3c81-4d26-aa31-cb027b9692da
begin
plot(Normal(mean(trait.(d25_p1.agents)), std(trait.(d25_p1.agents))), color=:black, label="Haploid")
plot!(Normal(mean(trait.(d25_p2.agents)), std(trait.(d25_p2.agents))), color=:blue, label="Diploid")
plot!(Normal(mean(trait.(d25_p4.agents)), std(trait.(d25_p4.agents))), color=:red, label="Tetraploid")
xlabel!("Phenotype before selection")
end

# ╔═╡ d9f04b11-8b31-4ad0-8486-16bf34e0c187
begin
plot(Normal(mean(trait.(d50_p1.agents)), std(trait.(d50_p1.agents))), color=:black, label="Haploid")
plot!(Normal(mean(trait.(d50_p2.agents)), std(trait.(d50_p2.agents))), color=:blue, label="Diploid")
plot!(Normal(mean(trait.(d50_p4.agents)), std(trait.(d50_p4.agents))), color=:red, label="Tetraploid")
xlabel!("Phenotype before selection")
end

# ╔═╡ c42be1ca-66a0-4abe-90f2-e72a37900b0b
begin
plot(Normal(mean(trait.(d75_p1.agents)), std(trait.(d75_p1.agents))), color=:black, label="Haploid")
plot!(Normal(mean(trait.(d75_p2.agents)), std(trait.(d75_p2.agents))), color=:blue, label="Diploid")
plot!(Normal(mean(trait.(d75_p4.agents)), std(trait.(d75_p4.agents))), color=:red, label="Tetraploid")
xlabel!("Phenotype before selection")
end

# ╔═╡ a9df72a6-b97b-4f77-9fde-faaeb123f9f3
trait.(d25_p1.agents)

# ╔═╡ Cell order:
# ╟─588b9cb0-a929-11eb-2ac2-e77ecad3b43f
# ╟─d2107ebd-72c0-48fb-8dbb-7a6478780104
# ╟─1a688141-c921-44c1-9685-818aff3c7b9e
# ╠═c1d263e3-d3d7-48f4-871f-05f82f858b8b
# ╠═8fa49926-3bf2-4eb9-95c8-634941600ea2
# ╠═c1c22e59-3fb6-45df-b793-3aa463398aec
# ╠═9e36a5a1-6461-4e8a-b088-cc838f86dd5e
# ╠═367df54e-3c81-4d26-aa31-cb027b9692da
# ╠═d9f04b11-8b31-4ad0-8486-16bf34e0c187
# ╠═c42be1ca-66a0-4abe-90f2-e72a37900b0b
# ╠═a9df72a6-b97b-4f77-9fde-faaeb123f9f3
# ╠═898ea668-2254-4927-8fcb-fdb1ccf3e6ec
# ╟─e7b41fc7-af1a-4b0e-aff7-f21bb85d1e9f
# ╠═a5f6483f-1567-4ebf-b58e-a08b229b2ef9
# ╟─e5596486-eab1-4a86-a31d-c33ded507add
# ╠═3098bc8d-42a5-4ff9-9931-9d7553a3b489
# ╠═ae903443-4310-45ce-beaa-53aec00db7ab
# ╠═5f227e3d-581c-4b50-aad9-92931427e3b2
# ╠═32f5be62-7fb2-400f-961e-8b394d91cb97
# ╠═e549e9e6-2d83-4fdb-967b-b668cc209687
# ╠═5b22b365-1898-4079-bd53-a164ed2f94fc
# ╟─9b6e587e-6266-42b8-853e-137d135422bd
# ╠═ded2e881-c5a2-4faa-83c0-0471e633a522
# ╠═7a021230-ccab-43bd-bef1-5ce25e9db559
# ╠═5b413f3d-f2dd-48fc-b761-1528b4604d7a
# ╟─1e6a7e16-3585-4135-a0fc-20e9c89105ea
# ╠═b1ccc66f-880d-40c5-8309-84e31c2fe1d1
# ╠═4bd0cf27-361a-4ec8-811e-f20f019c9c15
# ╠═81e0696e-1d12-410c-89e8-86e492810a62
