### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 2b840c70-8697-11eb-36e2-27dc0929a48b
using PolyStab

# ╔═╡ a74705e0-8699-11eb-1e70-77407cfbbcfb
using PolyStab: SimpleDeme, randagent, random_mating, allelefreqs

# ╔═╡ b79eb4e0-8696-11eb-0949-555c0c3c411f
md"""### Genetic variance in a mixed ploidy population"""

# ╔═╡ ded34c80-8699-11eb-33ce-214cb4f59699


# ╔═╡ a6337300-8699-11eb-1dfb-f7126f018d32
d = SimpleDeme(agents=randagent(0.2, 0.1, 100, 100))

# ╔═╡ Cell order:
# ╟─b79eb4e0-8696-11eb-0949-555c0c3c411f
# ╠═2b840c70-8697-11eb-36e2-27dc0929a48b
# ╠═a74705e0-8699-11eb-1e70-77407cfbbcfb
# ╠═ded34c80-8699-11eb-33ce-214cb4f59699
# ╠═a6337300-8699-11eb-1dfb-f7126f018d32
