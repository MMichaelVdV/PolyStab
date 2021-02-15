### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 49f46d76-0acf-11eb-32a2-3b0e0fbae9b6
using Random, Distributions, Plots, PlutoUI

# ╔═╡ 5ac23816-0acf-11eb-2a83-45fa592cecf1
md"""
# Single deme population genetics

I'll try a Pluto notebook, since I want to get acquainted with this new format (also I started to dislike the Jupyter notebooks). I'm trying to get down to the basics for a single deme dynamics here. This is in fact just an individual-based simulation of a haploid Wright-Fisher (WF) population.
"""

# ╔═╡ 8933cb92-0acf-11eb-31ff-37d2c3d86ecb
struct Agent{T}
    loci::Vector{T}
end

# ╔═╡ 87739a8e-0ad0-11eb-057b-0775eaa5aba8
begin
	randagent(d::Distribution, n) = Agent(rand(d, n))
	randagent(d::Distribution, n, N) = [Agent(rand(d, n)) for i=1:N]
end

# ╔═╡ 9be40b58-0acf-11eb-2577-654446a3c332
struct Deme{A,T}
    agents::Vector{A}
    K::Int64   # I'll use K instead of `carrying_capacity`...
    θ::T       # ... and θ instead of `optimal_phenotype`
end

# ╔═╡ f757bb24-0ad6-11eb-2cb4-df9cab6db0f6
md"""
As usual in `julia`, extending some methods from `Base` for our custom types when relevant makes life easier and code more readable (but beware, when extensions do not make much sense, for instance defining `Base.length` for something that is not array-like or so, it can make code less readable...)
"""

# ╔═╡ 2e59d926-0ad1-11eb-1717-517f339cdf74
begin
	Base.length(d::Deme) = length(d.agents)
	Base.length(a::Agent) = length(a.loci)
	Base.show(io::IO, d::Deme{A,T}) where {A,T} = 
		write(io, "Deme{$A,$T}(N=$(length(d)))")
	Base.push!(d::Deme, a::Agent) = push!(d.agents, a)
	Base.push!(a::Agent{T}, x::T) where T = push!(a.loci, x)
	Base.getindex(a::Agent, i) = a.loci[i]
	Base.getindex(d::Deme, i) = d.agents[i]
	Base.getindex(d::Deme, i, j) = d.agents[i][j]
	Base.rand(d::Deme) = rand(d.agents)
	Base.rand(d::Deme, n::Integer) = rand(d.agents, n)
end

# ╔═╡ 213694fa-0ad1-11eb-3fc5-6d177274c525
d = Deme(randagent(DiscreteNonParametric([0,1],[0.5,0.5]), 100, 100), 100, 0.)

# ╔═╡ 26ed8348-0ad3-11eb-346d-7f30ae5ee8c6
rand(d)

# ╔═╡ 5f94c272-0ad0-11eb-1eff-b71070c6732b
begin
	emptycopy(a::Agent{L}) where L = Agent(L[])
	emptycopy(d::Deme{A}) where A = Deme(A[], d.K, d.θ)
end

# ╔═╡ fc7db654-0ad3-11eb-1bc1-0757c681a7e4
md"""
Instead of empty copy and pushing to the array, it seems more natural (and ikely more efficient) to initialize an array of the right size and type, and filling it up.
"""

# ╔═╡ 18454a36-0ad5-11eb-2413-8f299ad04329
md"""
The function `random_mating` you defined is incorrect, well incorrect in the sense that it does not do what we usually understand under random mating. If we assume the Wright-Fisher model of random mating finite populations, the offspring generation is formed by randomly drawing parents *with replacement* from the current generation. You are drawing random parents without replacement. Moreover, what you have coded (in combination with your `recombine` function) amounts to doing haploid sexual reproduction (merging two genomes, doing recombination and reductive division) and adding both products of this process to the population, which I guess is relevant for some organisms, but is not a general or very pertinent abstract model of reproduction.

```julia
function random_mating(d::Deme)
    new_deme = emptycopy(d)
    nodes = d.agents
    while length(nodes) > 1
        shuffle!(nodes)
        n1 = pop!(nodes)
        n2 = pop!(nodes)
        nr = recombine(n1,n2)
        push!(new_deme,nr[1])
        push!(new_deme,nr[2])
    end
    new_deme
end
```
The following `mate` function picks a single offspring individual from a mating between two agents, and the random mating function samples parents with replacement. The resulting allel frequency dynamics are those of the Wright-Fisher model (and so at a single locus for instance the standard coalescent will apply).
"""

# ╔═╡ 6286d216-0ad7-11eb-131d-75a9434d668c
function mate(a::Agent{T}, b::Agent{T}) where T
    newloci = Vector{T}(undef, length(a))
    for i in 1:length(a)
		newloci[i] = rand() < 0.5 ? a[i] : b[i] 
    end
	Agent(newloci)
end 

# ╔═╡ e3149d3a-0acf-11eb-2569-47a2a0cf76a0
function random_mating(d::Deme{A}) where {A}
    newdeme = Vector{A}(undef, length(d))
    for i=1:length(d)
		newdeme[i] = mate(rand(d, 2)...)
	end
	Deme(newdeme, d.K, d.θ)
end 

# ╔═╡ ed04cf38-0ad7-11eb-1a8d-27f2a0e5355b
md"""
Your function for the heterozygosity computed the observed heterozygosity for the first locus right? I'll write it for all loci. I've been practicing on my functional programming style lately, so don't worry about the `map` and `mapreduce`, you can do it with for loops and ordinary sums as well, but just for fun.
"""

# ╔═╡ bd32afda-0ad7-11eb-1e3c-d3595ad626d7
function allelef(d::Deme)
	f = j->mapreduce(i->d[i,j], +, 1:length(d))/length(d)
	frequencies = map(f, 1:length(d[1]))
end

# ╔═╡ d794c158-0ae1-11eb-076e-51b7ea6617da
heterozygosities(d::Deme, freqs=allelef(d)) = map(p->p*(1-p), freqs)

# ╔═╡ de9b7e0a-0ad8-11eb-2c08-cb3e3f91d53a
heterozygosities(d)

# ╔═╡ 42053e30-0ada-11eb-1d72-cb994a1e0854
md"""
Note that the formula for the heterozygosity derived from the WF model is for the *expected* heterozygosity i.e. $\mathbb{E}[H_o(t)]$ where $H_o(t)$ is the *observed* heterozygosity, it is not a prediction of a deterministic model, but an expected value for a probabilistic model! So I won't call it `det_het` ;).
"""

# ╔═╡ 3c92aa5a-0ada-11eb-3155-5979c2c1c514
expected_heterozygosity(H₀, t, N) = ((1.0-1.0/N)^t)*H₀

# ╔═╡ 1aae55e6-0ae2-11eb-1fe9-53f9c6bce241
md"""
Now besides heterozygosities, we could also follow the additive genetic variance. Consider a haploid population with $n$ loci, if we assume each `1` allele has an additive effect of $a$ and the `0` allele an additive effect of 0, the expected contribution of locus $i$ to the phenotypic value $z$ is 

$E[X_i] = 0\times (1-p_i) + a\times p_i = a p_i$

Where we assume $p_i$ is the frequency of `1` at locus $i$. The expected phenotypic value $z$ is then $E[z] = \sum_i^n E[X_i]$ (note this is always true, also for linked loci). Recall that the variance is defined as

$Var(z) = E[(z - E[z])^2] = E[z^2 - 2E[z]z + E[z]^2]$
$\ = E[z^2] - 2E[z]E[z] + E[z]^2 = E[z^2] - E[z]^2$

and we have

$E[X_i^2] = 0 \times (1-p_i) + a^2 p_i = a^2  p_i$

so as a result 

$Var(X_i) = a^2 p_i - (a p_i)^2 = a^2 p_i(1-p_i)$

Assuming unlinked loci (i.e. probabilistic independence of the $X_i$), the additive genetic variance is then $Var(z) = \sum_i^n Var(X_i) = a^2 \sum_i^n p_i(1-p_i)$. Clearly, and expectedly, this is closely related to the heterozygosity. In fact the additive genetic variance in the unlinked case is a simple linear function of the average heterozygosity, so it is actually not interesting to follow seperately in these simulations. Note that this is straightforwardly extended to the case with locus-specific additive effects.
"""

# ╔═╡ 5e0a8c1a-0ad9-11eb-1d6d-bd3c5a983f18
function neutral_evolving_deme(d::Deme, ngen; fun=heterozygosities)
	stats = [fun(d)]
	for n=1:ngen
		d = random_mating(d)
		push!(stats, fun(d))
	end
	(stats=stats, deme=d, ngen=ngen)
end

# ╔═╡ b1aa2488-0adc-11eb-2da0-4d33dc9becff
md"""
- N $(@bind N Slider(1:100, default=50, show_value=true))
- L $(@bind L Slider(1:1000, default=100, show_value=true))
- p $(@bind p Slider(0.05:0.05:0.95, default=0.5, show_value=true))
- t $(@bind t Slider(1:1000, default=500, show_value=true))
"""

# ╔═╡ 6139af9a-0add-11eb-2baa-b1ca76a911e9
(N, L, p, t)

# ╔═╡ ee04a2b2-0ad9-11eb-00ce-47541d141b95
deme = Deme(randagent(DiscreteNonParametric([0,1],[p, 1-p]), L, N), N, 0.)

# ╔═╡ db9cb06a-0ad9-11eb-244b-693cf8cad0c5
sim = neutral_evolving_deme(deme, t)

# ╔═╡ 196ff97e-0ada-11eb-2336-85aa10e8975c
begin
	Hₒ = map(mean, sim.stats)
	plot(Hₒ, grid=false, color=:black, label="\$H_o(t)\$")
	plot!(0:sim.ngen+1, 
		t->expected_heterozygosity(p*(1-p), t, N),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 1cb76286-0ada-11eb-2c62-2152ffb909d2
md"""
Great, that's genetic drift as we'd expect it!
"""

# ╔═╡ Cell order:
# ╠═49f46d76-0acf-11eb-32a2-3b0e0fbae9b6
# ╟─5ac23816-0acf-11eb-2a83-45fa592cecf1
# ╠═8933cb92-0acf-11eb-31ff-37d2c3d86ecb
# ╠═87739a8e-0ad0-11eb-057b-0775eaa5aba8
# ╠═9be40b58-0acf-11eb-2577-654446a3c332
# ╟─f757bb24-0ad6-11eb-2cb4-df9cab6db0f6
# ╠═2e59d926-0ad1-11eb-1717-517f339cdf74
# ╠═213694fa-0ad1-11eb-3fc5-6d177274c525
# ╠═26ed8348-0ad3-11eb-346d-7f30ae5ee8c6
# ╠═5f94c272-0ad0-11eb-1eff-b71070c6732b
# ╠═fc7db654-0ad3-11eb-1bc1-0757c681a7e4
# ╟─18454a36-0ad5-11eb-2413-8f299ad04329
# ╠═6286d216-0ad7-11eb-131d-75a9434d668c
# ╠═e3149d3a-0acf-11eb-2569-47a2a0cf76a0
# ╟─ed04cf38-0ad7-11eb-1a8d-27f2a0e5355b
# ╠═bd32afda-0ad7-11eb-1e3c-d3595ad626d7
# ╠═d794c158-0ae1-11eb-076e-51b7ea6617da
# ╠═de9b7e0a-0ad8-11eb-2c08-cb3e3f91d53a
# ╟─42053e30-0ada-11eb-1d72-cb994a1e0854
# ╠═3c92aa5a-0ada-11eb-3155-5979c2c1c514
# ╟─1aae55e6-0ae2-11eb-1fe9-53f9c6bce241
# ╠═5e0a8c1a-0ad9-11eb-1d6d-bd3c5a983f18
# ╠═b1aa2488-0adc-11eb-2da0-4d33dc9becff
# ╠═6139af9a-0add-11eb-2baa-b1ca76a911e9
# ╠═ee04a2b2-0ad9-11eb-00ce-47541d141b95
# ╠═db9cb06a-0ad9-11eb-244b-693cf8cad0c5
# ╠═196ff97e-0ada-11eb-2336-85aa10e8975c
# ╟─1cb76286-0ada-11eb-2c62-2152ffb909d2
