### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 8ea08d1e-8bff-11eb-3fbb-f74cf53a2eed
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, DataFrames, StatsPlots

# ╔═╡ ded34c80-8699-11eb-33ce-214cb4f59699
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, evolving_neutraldeme, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy, trait_mean, randagent

# ╔═╡ b79eb4e0-8696-11eb-0949-555c0c3c411f
md"""### Genetic variance in a mixed ploidy population"""

# ╔═╡ 09e48960-9875-11eb-35fe-5122ae8dbe5f
md""" ##### Initialization of demes with different ploidy levels (haploid, diploid, tetraploid)"""

# ╔═╡ 217c5170-8e52-11eb-175c-c9b8476ae81b
d_p1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [1., 0., 0., 0.],100), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20.)

# ╔═╡ 3f70dcc0-8bfd-11eb-3df6-6709900e71ce
d_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [0., 1., 0., 0.],100), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20.)

# ╔═╡ 39b5e750-8dbd-11eb-3f64-c12fa3fe01d5
d_p4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [0., 0., 0., 1.],100), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20.)

# ╔═╡ 3106b1d0-9875-11eb-36ac-0d2b98a73fd4
md""" ##### Dynamics in a neutral deme (random mating, no selection)"""

# ╔═╡ a8d89200-8e52-11eb-1277-b1857f7ea0fe
begin
neutral_p1 = evolving_haploiddeme(d_p1, 500)
neutral_p2 = evolving_neutraldeme(d_p2, 500)
neutral_p4 = evolving_neutraldeme(d_p4, 500)
end

# ╔═╡ 66a25c90-8bff-11eb-2a53-09b98728f17c
begin
	Hₒ_p1 = map(mean, neutral_p1.het)
	Hₒ_p2 = map(mean, neutral_p2.het)
	Hₒ_p4 = map(mean, neutral_p4.het)
	expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	eh_1(H₀, t, N) = ((1.0-1.0/(N))^t)*H₀
	eh_2(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	eh_4(H₀, t, N) = ((1.0-1.0/(4*N))^t)*H₀
	plot(Hₒ_p1, grid=false, color=:black, label="\$H_op1(t)\$", title="LOH")
	plot!(Hₒ_p2, grid=false, color=:blue, label="\$H_op2(t)\$")
	plot!(Hₒ_p4, grid=false, color=:red, label="\$H_op4(t)\$")
	#plot!(0:sim_p2.ngen+1, 
		#t->expected_heterozygosity(0.5*(1-0.5), t, 50),
		#linestyle=:dash, color=:black, 
		#label = "\$(1-1/N)^t H_o(0)\$")
	plot!(0:neutral_p1.ngen+1, 
		t->eh_1(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	plot!(0:neutral_p2.ngen+1, 
		t->eh_2(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:blue, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	plot!(0:neutral_p4.ngen+1, 
		t->eh_4(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:red, 
		label = "\$(1-1/4N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 621f41a0-9885-11eb-2bf8-a530b5cc7e86
begin
n = 200 #number of loci
k = 2 #ploidy
#α = 1/√(k*n) 
α = .2 
a2 = randagent(0.5, α, n, 2)
b2 = randagent(0.5, α, n, 2)
o2 = map(x->trait(mate_p(a2,b2)), 1:2500) 

O_mean2 = mean(o2)
E_mean2 = 0.5*(trait(a2) + trait(b2)) 
O_var2 = var(o2)
E_var2 = n * (α^2 / 4)
	
a4 = randagent(0.5, α, n, 4)
b4 = randagent(0.5, α, n, 4)
o4 = map(x->trait(mate_p(a4,b4)), 1:5000) 

O_mean4 = mean(o4)
E_mean4 = 0.5*(trait(a4) + trait(b4)) 
O_var4 = var(o4)
E_var4 = 
	
plot(Normal(O_mean2, O_var2),color=:blue, label="Observed 2N")
plot!(Normal(O_mean4, O_var4),color=:red, label="Observed 4N")
plot!(Normal(E_mean2,E_var2),color=:black, label="Expected 2N")
xlabel!("\$phenotype\$")
ylabel!("\$\$")
end

# ╔═╡ 1b1541de-8dbf-11eb-041c-815c116efeeb
begin
	pop_p1 = map(mean, neutral_p1.p1)
	pop_p2 = map(mean, neutral_p2.p2)
	pop_p4 = map(mean, neutral_p4.p4)
	plot(pop_p1, grid=false, color=:black, label="\$pop_p1(t)\$", title="Population size")
	plot!(pop_p2, grid=false, color=:blue, label="\$pop_p2(t)\$")
	plot!(pop_p4, grid=false, color=:red, label="\$pop_p4(t)\$")

	xlabel!("\$t\$")
	ylims!(95,105)
	ylabel!("\$pop(t)\$")
end

# ╔═╡ 3fc09382-986b-11eb-29d6-57eeb0eddb77
md""" ##### Dynamics in a deme with densitiy dependence and stabilizing selection"""

# ╔═╡ 47dc60c0-9876-11eb-3d70-cd03509cbad8
begin
stabsel_p1 = evolving_selectiondeme(d_p1,500)
stabsel_p2 = evolving_selectiondeme(d_p2,500)
stabsel_p4 = evolving_selectiondeme(d_p4,500)
end

# ╔═╡ c1798110-9885-11eb-2e75-ef20122e8add
md""" Stabilizing selection is expected to further degrade genetic variance."""

# ╔═╡ fc96ac10-9884-11eb-1735-d751fd99c805
begin
	Hₒs_p1 = map(mean, stabsel_p1.het)
	Hₒs_p2 = map(mean, stabsel_p2.het)
	Hₒs_p4 = map(mean, stabsel_p4.het)

	plot(Hₒs_p1, grid=false, color=:black, label="\$H_op1(t)\$", title="LOH")
	plot!(Hₒs_p2, grid=false, color=:blue, label="\$H_op2(t)\$")
	plot!(Hₒs_p4, grid=false, color=:red, label="\$H_op4(t)\$")
	#plot!(0:sim_p2.ngen+1, 
		#t->expected_heterozygosity(0.5*(1-0.5), t, 50),
		#linestyle=:dash, color=:black, 
		#label = "\$(1-1/N)^t H_o(0)\$")
	plot!(0:neutral_p1.ngen+1, 
		t->eh_1(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	plot!(0:neutral_p2.ngen+1, 
		t->eh_2(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:blue, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	plot!(0:neutral_p4.ngen+1, 
		t->eh_4(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:red, 
		label = "\$(1-1/4N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 66257e10-9879-11eb-2367-d706689f326a
begin
	traitmean_p1 = map(mean, stabsel_p1.tm)
	traitmean_p2 = map(mean, stabsel_p2.tm)
	traitmean_p4 = map(mean, stabsel_p4.tm)
	plot(traitmean_p1, grid=false, color=:black, label=false)
	plot!(traitmean_p2, grid=false, color=:blue, label=false)
	plot!(traitmean_p4, grid=false, color=:red, label=false)
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 63ad1c50-987f-11eb-1512-df6fc858aaca
begin
traitp1 = plot(traitmean_p1, grid=false, color=:black, label=false)
traitp2 = plot(traitmean_p2, grid=false, color=:blue, label=false)
traitp3 = plot(traitmean_p4, grid=false, color=:red, label=false)
plot(traitp1,traitp2,traitp3, layout = (3, 1))
end

# ╔═╡ c9b2b8d0-9879-11eb-294a-952184f8bf9c
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

# ╔═╡ 9cdac49e-9884-11eb-2be7-832739ca5d66
begin
	popp1 = plot(popsize_p1, grid=false, color=:black, label=false)
	popp2 = plot(popsize_p2, grid=false, color=:blue, label=false)
	popp4 = plot(popsize_p4, grid=false, color=:red, label=false)
	plot(popp1,popp2,popp4, layout = (3, 1))
end

# ╔═╡ 21c69f40-9876-11eb-0f98-cf5c73410b72
md""" ###### Effect of genotype to phenotype map"""

# ╔═╡ 487c7150-9876-11eb-0bf3-4557f19a59af


# ╔═╡ Cell order:
# ╟─b79eb4e0-8696-11eb-0949-555c0c3c411f
# ╠═8ea08d1e-8bff-11eb-3fbb-f74cf53a2eed
# ╠═ded34c80-8699-11eb-33ce-214cb4f59699
# ╟─09e48960-9875-11eb-35fe-5122ae8dbe5f
# ╠═217c5170-8e52-11eb-175c-c9b8476ae81b
# ╠═3f70dcc0-8bfd-11eb-3df6-6709900e71ce
# ╠═39b5e750-8dbd-11eb-3f64-c12fa3fe01d5
# ╟─3106b1d0-9875-11eb-36ac-0d2b98a73fd4
# ╠═a8d89200-8e52-11eb-1277-b1857f7ea0fe
# ╠═66a25c90-8bff-11eb-2a53-09b98728f17c
# ╠═621f41a0-9885-11eb-2bf8-a530b5cc7e86
# ╠═1b1541de-8dbf-11eb-041c-815c116efeeb
# ╟─3fc09382-986b-11eb-29d6-57eeb0eddb77
# ╠═47dc60c0-9876-11eb-3d70-cd03509cbad8
# ╟─c1798110-9885-11eb-2e75-ef20122e8add
# ╟─fc96ac10-9884-11eb-1735-d751fd99c805
# ╟─66257e10-9879-11eb-2367-d706689f326a
# ╟─63ad1c50-987f-11eb-1512-df6fc858aaca
# ╟─c9b2b8d0-9879-11eb-294a-952184f8bf9c
# ╟─9cdac49e-9884-11eb-2be7-832739ca5d66
# ╟─21c69f40-9876-11eb-0f98-cf5c73410b72
# ╠═487c7150-9876-11eb-0bf3-4557f19a59af
