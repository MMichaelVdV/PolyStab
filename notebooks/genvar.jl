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
d_p = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 1., 200, [2], 50, d=2.),randagent_p(0.5, 0.1, 250, [4], 0, d=4.)), OV = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.] )

# ╔═╡ 6cea6720-8bfd-11eb-3640-d5d0098ebd7a
sim = neutral_evolving_deme(d_p, 500)

# ╔═╡ 66a25c90-8bff-11eb-2a53-09b98728f17c
begin
	Hₒ = map(mean, sim.het)
	expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	plot(Hₒ, grid=false, color=:black, label="\$H_o(t)\$")
	plot!(0:sim.ngen+1, 
		t->expected_heterozygosity(0.5*(1-0.5), t, 50),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ c9e15e00-8c04-11eb-007d-ebe495b92603
recombine_poly(d_p)

# ╔═╡ a68923f2-8c06-11eb-131d-133c229a5505
begin
	dp = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 1., 200, [2], 100, d=2.),randagent_p(0.5, 0.1, 250, [4], 0, d=4.)), OV = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.] )
	dps = foldl((dp,i)->random_mating_mixedp(dp), 1:1000, init=dp)
	sort(unique(allelefreqs_p(dps))) == [0., 1.]
end

# ╔═╡ Cell order:
# ╟─b79eb4e0-8696-11eb-0949-555c0c3c411f
# ╠═8ea08d1e-8bff-11eb-3fbb-f74cf53a2eed
# ╠═2b840c70-8697-11eb-36e2-27dc0929a48b
# ╠═ded34c80-8699-11eb-33ce-214cb4f59699
# ╠═3f70dcc0-8bfd-11eb-3df6-6709900e71ce
# ╠═6cea6720-8bfd-11eb-3640-d5d0098ebd7a
# ╠═66a25c90-8bff-11eb-2a53-09b98728f17c
# ╠═c9e15e00-8c04-11eb-007d-ebe495b92603
# ╠═a68923f2-8c06-11eb-131d-133c229a5505
