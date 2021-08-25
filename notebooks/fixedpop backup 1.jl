### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ b107f00a-8c3a-446a-a1d0-32a1323cfb6a
using StatsBase, Plots

# ╔═╡ 07704d0b-0cbc-4b44-8fed-1cd8831ad075
using PolyStab: AbstractDeme, randagent_p, MixedPloidyDeme, malthusian_fitness, number_of_offspring, ismock, mate_p, ploidy_freq

# ╔═╡ 8d201245-5a5d-44ef-a7af-e04dc92f2f30
using PolyStab: heterozygosities_p, trait_mean, allelefreqs_p, f_trait_agents, mutate

# ╔═╡ e6a288a2-b8f2-11eb-2d78-2d5dc5c6684b
md""" ##### Some implementations of selection methods with constant population size."""

# ╔═╡ ced90b8a-3f8a-4b22-8ff8-f143e312d782
md""" Population size in a random mating deme without any form of selection crashes in mixed ploidy system because of cytoload"""

# ╔═╡ 7b8f7657-8b8a-4dbe-8ec0-746a67555450
begin

"""
	mating_PnB_x(d::MixedPloidyDeme{A})

Mating in a mixed ploidy deme with unreduced gamete formation and partner
choice weighted by malthusian fitness.
"""
function mating_PnB(d::AbstractDeme{A}) where A
	new_agents = A[]
	fitnesses = exp.(malthusian_fitness(d))
	while length(new_agents) < d.K
		B1 = sample(d.agents)
		noff = number_of_offspring(d, B1)
        Bs = sample(d.agents, weights(fitnesses), noff) 
        offspring = filter(!ismock, map(B2->mate_p(B1, B2, d), Bs))
        new_agents = vcat(new_agents, offspring)
    end
    d(new_agents)
end
	
	
	"""
		function truncation(fitness, N::Int)
	
	Return top `N` of highest fitness inividuals.
	"""
	function truncation(fitness, N)
    λ = length(fitness)
    @assert λ >= N "Cannot select more then $(λ) elements"
    idx = reverse(sortperm(fitness))
    return idx[1:N]
	end
	
	
"""
	mating_PnB(d::AbstractDeme{A})
	
Mating in a mixed ploidy deme.
"""
	function mating_PnB(d::AbstractDeme{A}, T) where A
		top_agents =  A[]
		new_agents =  A[]
		fitnesses = exp.(malthusian_fitness(d))
		
		N = round(Int64, T*length(d))
		
		for i in truncation(fitnesses, N)
			push!(top_agents, d[i])
		end
		deme_top = d(top_agents)
		
		fitnesses_top = exp.(malthusian_fitness(deme_top))
		
		while length(new_agents) < d.K
	 		B1 = sample(top_agents)
			noff = number_of_offspring(d,B1) 
			B2 = sample(top_agents, weights(fitnesses_top))
			m = mate_p(B1,B2,d)
			if m != 0
				for c in 1:noff 
					push!(new_agents, m)
				end
			end
		end
	    d(new_agents)
	end
	
end

# ╔═╡ 7f65aed1-cf90-44f9-b6db-766e5ff26410
begin


"""
	evolving_selectiondeme(d::AbstractDeme, ngen)

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.
"""
function evolving_selectiondeme(d::MixedPloidyDeme, ngen; 
	heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = 
	allelefreqs_p, pf = ploidy_freq, fta = f_trait_agents)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[2]]
	p3 = [ploidy_freq(d)[3]]
	p4 = [ploidy_freq(d)[4]]
	fta = [f_trait_agents(d)]
	fit = [malthusian_fitness(d)]
	
	for n=1:ngen
		d = mating_PnB(d)
		d = mutate(d) 
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[2])
		push!(p3, ploidy_freq(d)[3])
		push!(p4, ploidy_freq(d)[4])
		push!(fta, f_trait_agents(d))
		push!(fit, malthusian_fitness(d))
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, het=het,tm=tm, af=af, fta=fta, fit=fit) 
end

"""
	evolving_selectiondeme(d::AbstractDeme, ngen)

Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
gamete formation.
"""
function evolving_selectiondeme(d::MixedPloidyDeme, T, ngen; 
	heterozygosities_p=heterozygosities_p, fit=malthusian_fitness, trait_mean = trait_mean, allelefreqs_p = 
	allelefreqs_p, pf = ploidy_freq, fta = f_trait_agents)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[2]]
	p3 = [ploidy_freq(d)[3]]
	p4 = [ploidy_freq(d)[4]]
	fta = [f_trait_agents(d)]
	fit = [malthusian_fitness(d)]
	
	for n=1:ngen
		d = mating_PnB(d, T)
		d = mutate(d) 
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[2])
		push!(p3, ploidy_freq(d)[3])
		push!(p4, ploidy_freq(d)[4])
		push!(fta, f_trait_agents(d))
		push!(fit, malthusian_fitness(d))
		
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, het=het,tm=tm, af=af, fta=fta, fit=fit) 
end

end

# ╔═╡ 19950de2-7e55-4081-b1d2-c40738958079
d_mp = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],200), UG = [0.0 0.0 0.0 0.0; 1. 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], K=200)

# ╔═╡ 45d3e279-8261-47b5-9928-b6b284df44da
malthusian_fitness(d_mp)

# ╔═╡ 7659cd05-db92-4c59-a867-b3dc2e2849b4
ss = foldl((d_mp,i)->mating_PnB(d_mp), 1:500, init=d_mp)

# ╔═╡ b623b469-f767-4f1b-a326-86c0de8055c4
mean(malthusian_fitness(ss))

# ╔═╡ 7f2a45b6-32c4-4d3c-9291-f0d4b392c8ba
ploidy_freq(ss)

# ╔═╡ 82608dc5-b936-495f-a142-b107d1b48500
ts = foldl((d_mp,i)->mating_PnB(d_mp, 0.05), 1:500, init=d_mp)

# ╔═╡ a857eca1-e42c-4f37-bd90-9a1306c07c37
mean(malthusian_fitness(ts))

# ╔═╡ bb30fc43-0222-4899-ba33-6f78f9242058
ploidy_freq(ts)

# ╔═╡ 0b11f4fd-36a0-4f01-911a-b9c2ae403899
begin
	ss_p2 = evolving_selectiondeme(d_mp, 250)
	ts_p2 = evolving_selectiondeme(d_mp, 0.25, 250)
end

# ╔═╡ bc1497b5-1a4f-41de-b4ad-06a0fd3da6ff
ss_p2

# ╔═╡ 7b109a30-85c0-4e22-8f32-49ff0bae0cb7
begin
	traitmean_ssp2 = map(mean, ss_p2.tm)
	
	stabsel = plot(traitmean_ssp2, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright, title="Stabsel")
	
	for (i,t) in enumerate(ss_p2.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_mp.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	#ylims!(18,22)
end

# ╔═╡ b7e56a22-acba-41a4-bee2-0dce1df84469
begin
	
	traitmean_tsp2 = map(mean, ts_p2.tm)
	
	truncsel = plot(traitmean_tsp2, grid=false, color=:green, label=false,linewidth=3,legend=:bottomright, title="Truncsel")
	
	for (i,t) in enumerate(ts_p2.fta)
	scatter!([i for x in 1:10],t,label=false,colour="green",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_mp.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	#ylims!(18,22)
end

# ╔═╡ 2d0a84af-d163-489f-bf26-1df9d173f67d
begin
	fitmean_ssp2 = map(mean, ss_p2.fit)
	
	plot(fitmean_ssp2, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright, title="Stabsel")
	
	for (i,t) in enumerate((ss_p2.fit))
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_mp.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(-10,1)
end

# ╔═╡ 29aa3088-99be-432e-8213-7301297eac16
begin
	fitmean_tsp2 = map(mean, ts_p2.fit)
	
	plot(fitmean_tsp2, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright, title="Truncsel")
	
	for (i,t) in enumerate((ts_p2.fit))
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_mp.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(-10,1)
end

# ╔═╡ 5a6d1578-1321-49b7-a9a8-306206ac1824
begin
	plot(ss_p2.pop, grid=false, color=:black, label=false)
	plot!(ss_p2.p2, grid=false, color=:red, label=false)
	plot!(ss_p2.p4, grid=false, color=:blue, label=false)
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ dfb850cb-6b61-4143-b60b-fada07cc02a7
begin
	plot(ts_p2.pop, grid=false, color=:black, label=false)
	plot!(ts_p2.p2, grid=false, color=:red, label=false)
	plot!(ts_p2.p4, grid=false, color=:blue, label=false)
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ b6d487e5-a141-4244-a98e-36b99335afaf
plot(stabsel, truncsel)

# ╔═╡ Cell order:
# ╟─e6a288a2-b8f2-11eb-2d78-2d5dc5c6684b
# ╟─ced90b8a-3f8a-4b22-8ff8-f143e312d782
# ╠═b107f00a-8c3a-446a-a1d0-32a1323cfb6a
# ╠═07704d0b-0cbc-4b44-8fed-1cd8831ad075
# ╠═7b8f7657-8b8a-4dbe-8ec0-746a67555450
# ╠═7f65aed1-cf90-44f9-b6db-766e5ff26410
# ╠═19950de2-7e55-4081-b1d2-c40738958079
# ╠═45d3e279-8261-47b5-9928-b6b284df44da
# ╠═7659cd05-db92-4c59-a867-b3dc2e2849b4
# ╠═b623b469-f767-4f1b-a326-86c0de8055c4
# ╠═7f2a45b6-32c4-4d3c-9291-f0d4b392c8ba
# ╠═82608dc5-b936-495f-a142-b107d1b48500
# ╠═a857eca1-e42c-4f37-bd90-9a1306c07c37
# ╠═bb30fc43-0222-4899-ba33-6f78f9242058
# ╠═8d201245-5a5d-44ef-a7af-e04dc92f2f30
# ╠═0b11f4fd-36a0-4f01-911a-b9c2ae403899
# ╠═bc1497b5-1a4f-41de-b4ad-06a0fd3da6ff
# ╠═7b109a30-85c0-4e22-8f32-49ff0bae0cb7
# ╠═b7e56a22-acba-41a4-bee2-0dce1df84469
# ╠═2d0a84af-d163-489f-bf26-1df9d173f67d
# ╠═29aa3088-99be-432e-8213-7301297eac16
# ╠═5a6d1578-1321-49b7-a9a8-306206ac1824
# ╠═dfb850cb-6b61-4143-b60b-fada07cc02a7
# ╠═b6d487e5-a141-4244-a98e-36b99335afaf
