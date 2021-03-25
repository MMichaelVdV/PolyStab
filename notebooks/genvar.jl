### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 8ea08d1e-8bff-11eb-3fbb-f74cf53a2eed
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes

# ╔═╡ 2b840c70-8697-11eb-36e2-27dc0929a48b
using PolyStab

# ╔═╡ ded34c80-8699-11eb-33ce-214cb4f59699
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_deme_UG, evolving_deme_ploidyvar, heterozygosities_p, allelefreqs_p, neutral_evolving_deme, recombine_poly, random_mating_mixedp

# ╔═╡ b79eb4e0-8696-11eb-0949-555c0c3c411f
md"""### Genetic variance in a mixed ploidy population"""

# ╔═╡ 3f70dcc0-8bfd-11eb-3df6-6709900e71ce
d_p2 = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 1., 200, [2], 500, d=2.),randagent_p(0.5, 0.1, 250, [4], 0, d=4.)), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.] )

# ╔═╡ 6cea6720-8bfd-11eb-3640-d5d0098ebd7a
sim_p2 = neutral_evolving_deme(d_p2, 500)

# ╔═╡ 39b5e750-8dbd-11eb-3f64-c12fa3fe01d5
d_p4 = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 1., 200, [2], 0, d=2.),randagent_p(0.5, 0.1, 250, [4], 500, d=4.)), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.] )

# ╔═╡ 4dddd0d0-8dbd-11eb-37ae-d3b33fb87952
sim_p4 = neutral_evolving_deme(d_p4, 500)

# ╔═╡ bc1ee9d0-8dbd-11eb-2d41-af9fd9e2cce4
heterozygosities_p(d_p2)

# ╔═╡ c9616230-8dbd-11eb-27df-27dbfd045c89
heterozygosities_p(d_p4)

# ╔═╡ 66a25c90-8bff-11eb-2a53-09b98728f17c
begin
	Hₒ_p2 = map(mean, sim_p2.het)
	Hₒ_p4 = map(mean, sim_p4.het)
	expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	plot(Hₒ_p2, grid=false, color=:blue, label="\$H_op2(t)\$")
	plot!(Hₒ_p4, grid=false, color=:red, label="\$H_op4(t)\$")
	plot!(0:sim_p2.ngen+1, 
		t->expected_heterozygosity(0.5*(1-0.5), t, 500),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ a68923f2-8c06-11eb-131d-133c229a5505
begin
dp = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 1., 200, [2], 100, d=2.),randagent_p(0.5, 0.1, 250, [4], 0, d=4.)), OV = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.] )
dps = foldl((dp,i)->random_mating_mixedp(dp), 1:1000, init=dp)
sort(unique(allelefreqs_p(dps))) == [0., 1.]
end

# ╔═╡ 1b1541de-8dbf-11eb-041c-815c116efeeb
begin
	pop_p2 = map(mean, sim_p2.p2)
	pop_p4 = map(mean, sim_p4.p4)
	plot(pop_p2, grid=false, color=:blue, label="\$pop_p2(t)\$")
	plot!(pop_p4, grid=false, color=:red, label="\$pop_p4(t)\$")

	xlabel!("\$t\$")
	ylabel!("\$pop(t)\$")
end

# ╔═╡ Cell order:
# ╟─b79eb4e0-8696-11eb-0949-555c0c3c411f
# ╠═8ea08d1e-8bff-11eb-3fbb-f74cf53a2eed
# ╠═2b840c70-8697-11eb-36e2-27dc0929a48b
# ╠═ded34c80-8699-11eb-33ce-214cb4f59699
# ╠═3f70dcc0-8bfd-11eb-3df6-6709900e71ce
# ╠═6cea6720-8bfd-11eb-3640-d5d0098ebd7a
# ╠═39b5e750-8dbd-11eb-3f64-c12fa3fe01d5
# ╠═4dddd0d0-8dbd-11eb-37ae-d3b33fb87952
# ╠═bc1ee9d0-8dbd-11eb-2d41-af9fd9e2cce4
# ╠═c9616230-8dbd-11eb-27df-27dbfd045c89
# ╠═66a25c90-8bff-11eb-2a53-09b98728f17c
# ╠═a68923f2-8c06-11eb-131d-133c229a5505
# ╠═1b1541de-8dbf-11eb-041c-815c116efeeb
