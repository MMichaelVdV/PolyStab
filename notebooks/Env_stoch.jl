### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 902601f0-8c7b-4070-b628-101df46d5346
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, DataFrames, StatsPlots

# ╔═╡ cf436571-d2d2-4430-bb97-3886d3216993
using PolyStab: ploidy_freq, f_trait_agents, mating_PnB, mutate

# ╔═╡ 6d8f35f2-56aa-47c1-803d-244595383f18
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, evolving_neutraldeme, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy, trait_mean, randagent, evolving_haploiddeme, mate_p, evolving_selectiondemeh

# ╔═╡ 1bb8e789-990c-49fc-9d20-34cb106b32cc
using PolyStab: malthusian_fitness, trait_add, trait_dom, trait_rec

# ╔═╡ e3bd0be8-3fa2-433f-82a4-bfc8da8332b1
using DifferentialEquations

# ╔═╡ 1cc95f70-d78a-11eb-2d30-39281eef250d
md""" #### Environmental stochasticity """

# ╔═╡ f0b619db-547f-4f5d-a158-4a2ab8d19402


# ╔═╡ b6f72f70-8069-4ade-a94f-479e85df444b
"""
    trait_exp(a::Agent)
"""
trait_exp(a::Agent) = sum(a)


# ╔═╡ 58bb47de-128a-445d-8293-84b2d77bc045
"""
	evolving_selectiondeme(d::AbstractDeme, ngen)

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.
"""
function evolving_selectiondemeenv(d::MixedPloidyDeme, ftgmap, env, ngen; 
heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = 
allelefreqs_p, pf = ploidy_freq, fta = f_trait_agents)
het = [heterozygosities_p(d)]
pop = [length(d)]
tm = [trait_mean(d, ftgmap)]
af = [allelefreqs_p(d)]
p2 = [ploidy_freq(d)[2]]
p3 = [ploidy_freq(d)[3]]
p4 = [ploidy_freq(d)[4]]
fta = [f_trait_agents(d, ftgmap)]

for n=1:ngen
	d = MixedPloidyDeme(agents=d.agents, K=d.K, θ=env[n], rm=d.rm, Vs=d.Vs, α=d.α, μ=d.μ, OV=d.OV, UG=d.UG)
	d = mating_PnB_x(d, ftgmap)
	d = mutate(d) 
	push!(het, heterozygosities_p(d))
	push!(pop, length(d))
	push!(tm, trait_mean(d, ftgmap))
	push!(af, allelefreqs_p(d))
	push!(p2, ploidy_freq(d)[2])
	push!(p3, ploidy_freq(d)[3])
	push!(p4, ploidy_freq(d)[4])
	push!(fta, f_trait_agents(d, ftgmap))
	
end
(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, het=het,tm=tm, af=af, fta=fta) 
end

# ╔═╡ 3ac3447f-0fa3-4beb-8662-8635d822fd7f
envs = 0.01 .* sin.(1/50 .* (1:500)) .+ 20

# ╔═╡ fdd86050-c74c-461e-af13-3296500021c3
env = [20. for in in 1:500]

# ╔═╡ fb148a79-513d-44bf-b3c7-3783e7fd2a29
plot(env)

# ╔═╡ 1a6d7b53-7e51-483a-aca5-4616a86c3ba1
d_p2env = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [0., 0., 0., 1.],200), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 0.1, α = 0.4)

# ╔═╡ ef0b8f99-3f58-46eb-a090-4d23d33be5f1
trait.(d_p2env.agents)

# ╔═╡ 599d3e47-82a9-49f2-8825-b37c2882b6c6
stabselenv_p2 = evolving_selectiondemeenv(d_p2env, trait_exp, env, 500)

# ╔═╡ 78e6817d-1133-4ee5-b1bc-30b405351bcd
begin
traitmean_p2env = map(mean, stabselenv_p2.tm)
popsize_p2env = map(mean, stabselenv_p2.pop)
end

# ╔═╡ c033b4fd-c58f-4010-ad58-f9b728040d5c
begin
	p2env = plot(traitmean_p2env, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright, title="Diploid")
	for (i,t) in enumerate(stabselenv_p2.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_p2env.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(18,22)
end

# ╔═╡ 074d4e74-eaf8-4965-aff2-004fb705e650
begin
	plot(stabselenv_p2.pop, grid=false, color=:black, label=false)
	plot!(stabselenv_p2.p2, grid=false, color=:red, label=false)
	plot!(stabselenv_p2.p4, grid=false, color=:blue, label=false)
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ 93a1cd1c-a456-4fa1-9b41-65b4b4682c79
md""" ### Addendum: Diffusion equations"""

# ╔═╡ 2ccae953-4beb-419a-9cc3-c5a06c303a2c
begin
	α=1
	β=1
	u₀=1/2
	f(u,p,t) = α*u
	g(u,p,t) = β*u
	dt = 1//2^(4)
	tspan = (0.0,1.0)
	prob = SDEProblem(f,g,u₀,(0.0,1.0))
end

# ╔═╡ 1069ee6a-da36-4183-b93c-79d4e6b2e1d4
begin
	sol = [solve(prob,EM(),dt=dt) for i in 1:10]
	plotly() # Using the Plotly backend
	plot(sol)
end

# ╔═╡ Cell order:
# ╟─1cc95f70-d78a-11eb-2d30-39281eef250d
# ╠═902601f0-8c7b-4070-b628-101df46d5346
# ╠═cf436571-d2d2-4430-bb97-3886d3216993
# ╠═f0b619db-547f-4f5d-a158-4a2ab8d19402
# ╠═6d8f35f2-56aa-47c1-803d-244595383f18
# ╠═1bb8e789-990c-49fc-9d20-34cb106b32cc
# ╠═b6f72f70-8069-4ade-a94f-479e85df444b
# ╠═58bb47de-128a-445d-8293-84b2d77bc045
# ╠═3ac3447f-0fa3-4beb-8662-8635d822fd7f
# ╠═fdd86050-c74c-461e-af13-3296500021c3
# ╠═fb148a79-513d-44bf-b3c7-3783e7fd2a29
# ╠═1a6d7b53-7e51-483a-aca5-4616a86c3ba1
# ╠═ef0b8f99-3f58-46eb-a090-4d23d33be5f1
# ╠═599d3e47-82a9-49f2-8825-b37c2882b6c6
# ╠═78e6817d-1133-4ee5-b1bc-30b405351bcd
# ╠═c033b4fd-c58f-4010-ad58-f9b728040d5c
# ╠═074d4e74-eaf8-4965-aff2-004fb705e650
# ╠═93a1cd1c-a456-4fa1-9b41-65b4b4682c79
# ╠═e3bd0be8-3fa2-433f-82a4-bfc8da8332b1
# ╠═2ccae953-4beb-419a-9cc3-c5a06c303a2c
# ╠═1069ee6a-da36-4183-b93c-79d4e6b2e1d4
