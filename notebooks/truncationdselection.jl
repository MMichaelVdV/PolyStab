### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 748df908-da3d-4025-8a0c-369678f82bae
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, StatsPlots

# ╔═╡ 7b60ca0b-9ae9-4b41-bd6f-744372135a7a
using PolyStab: Agent, IslandDeme, randagent_p, heterozygosities_p, trait_mean, allelefreqs_p, ploidy_freq, f_trait_agents, mutate, trait, popsize, number_of_offspring, mate_p

# ╔═╡ c016ac10-b047-11eb-37a1-f198dc39687e
md""" ### Truncation selection"""

# ╔═╡ 985ddc8f-d065-449b-9eeb-41003b7449a8
md""" One model has received particular attention in the context of genetic loads. This is a  threshold, or truncation selection model. It is assumed that individuals with a  certain number of mutants or more are eliminated while those with less than this number are retained. With strict truncation at a  threshold of m mutants, then each eliminated individual has at least m mutants, and the considerations mentioned above become applicable. This kind of model has been discussed by Kimura and Maruyama (1966) and especially by King (1966) who empha-sized that the mutation load can be considerably reduced in such a situation. From Crow (1970). """

# ╔═╡ a1df4b9d-2ab3-4d20-9399-778330ed200c
md""" Truncation selection is described by either the fraction, p, of the population saved or the threshold phenotypic value, T, below (or above) which individuals are culled (Walsh and Lynch, 2018). The implementation here uses the fraction p of the population."""

# ╔═╡ d9f6864f-fdb6-43b8-b8bb-3ea98af682c7
begin
	"""
		directional_selection(d::AbstractDeme,a::Agent)
	"""
	function directional_selection(d::IslandDeme,a::Agent)
	    N = length(d)
	    z = trait(a)
	    return d.β*(z-d.θ)
	end 
	
	
	"""
		directional_selection(d::IslandDeme)
	"""
	function directional_selection(d::IslandDeme)
	    N = length(d)
	    fitnesses = Float64[]
		for agent in d.agents
			z = trait(agent)
	    	f = d.β*(z-d.θ)
			push!(fitnesses, f)
		end
		return fitnesses
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
		mating_PnB(d::IslandDeme{A})
	
	Mating in a mixed ploidy islanddeme with unreduced gamete formation directional selection.
	"""
	function mating_PnB(d::IslandDeme{A}, T) where A
		top_agents =  A[]
		new_agents =  A[]
		fitnesses = exp.(directional_selection(d))
		
		if length(d) < 10
			N = length(d)
		else N = round(Int64, T*length(d))
		end
		
		for i in truncation(fitnesses, N)
			push!(top_agents, d[i])
		end
		deme_top = d(top_agents)
		
		fitnesses_top = exp.(directional_selection(deme_top))
		
		for j in 1:length(top_agents)
	 		B1 = top_agents[j]
			noff = popsize(d) < 350 ? number_of_offspring(d,B1) : 1 #limits popsize to prevent exponential explosion
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
	
	
	"""
		evolving_islanddeme(d::AbstractDeme, ngen)
	
	Simulate a single deme with mixed ploidy, malthusian fitness and unreduced
	gamete formation.
	"""
	function evolving_islanddeme(d::IslandDeme, T, ngen; 
		heterozygosities_p=heterozygosities_p, fit=directional_selection, trait_mean = trait_mean, allelefreqs_p = 
		allelefreqs_p, pf = ploidy_freq, fta = f_trait_agents)
		het = [heterozygosities_p(d)]
		pop = [length(d)]
		tm = [trait_mean(d)]
		af = [allelefreqs_p(d)]
		p2 = [ploidy_freq(d)[2]]
		p3 = [ploidy_freq(d)[3]]
		p4 = [ploidy_freq(d)[4]]
		fta = [f_trait_agents(d)]
		
		for n=1:ngen
			d = mating_PnB(d, T)
			#d = mutate(d) 
			push!(het, heterozygosities_p(d))
			push!(pop, length(d))
			push!(tm, trait_mean(d))
			push!(af, allelefreqs_p(d))
			push!(p2, ploidy_freq(d)[2])
			push!(p3, ploidy_freq(d)[3])
			push!(p4, ploidy_freq(d)[4])
			push!(fta, f_trait_agents(d))
			
		end
		(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=ngen, het=het,tm=tm, af=af, fta=fta) 
	end
end

# ╔═╡ 11a96686-826f-4875-aec4-2f7d619bc133
islanddeme = IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], β=0.25, θ=15.5)

# ╔═╡ 15902758-b8bd-4b9d-b0d9-46e92007972c
migrant = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1)[1]

# ╔═╡ d98d3d38-2352-4800-89f4-b077b819d040
trait(migrant)

# ╔═╡ 2342c60c-7e5a-47a3-b23d-a8b8fca19d7c
push!(islanddeme.agents,migrant)

# ╔═╡ ded1736f-e18e-4174-9a6b-7df33b804c9e
begin
island_p2 = evolving_islanddeme(islanddeme, .85, 30)
while popsize(island_p2.deme) < 5
	island_p2 = evolving_islanddeme(islanddeme, 0.85, 30)
	end
end


# ╔═╡ 2eaa35ba-f42e-4417-82d9-7dc006070d8c
begin
	pf2_p2 = island_p2.p2
	pf3_p2 = island_p2.p3
	pf4_p2 = island_p2.p4
	p21 = plot(pf2_p2, grid=false, color=:blue, label="diploids", ylims = (1,250), legend=:bottomright)
	#plot!(pf3_p1, grid=false, color=:green, label="triploids")
	plot!(pf4_p2, grid=false, color=:red, label="tetraploids")
	hline!([islanddeme.K],label ="K",colour = "black",linestyle=:dash)
	
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ a7ce43d2-9ccd-4813-bcd5-b405eb5c0008
begin
	traitmean_ploidy2 = map(mean, island_p2.tm)
	p22 = plot(traitmean_ploidy2, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright) #title="Stabilizing selection"
	for (i,t) in enumerate(island_p2.fta)
	scatter!([i for x in 1:length(t)],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([islanddeme.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ 83cc65c3-3a5f-47fc-8055-7e319000f7da
begin
	Hₒ_ploidyhet2 = map(mean, 2 .*island_p2.het)
	expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	p23 = plot(Hₒ_ploidyhet2, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH", yscale = :log)
	plot!(1:island_p2.ngen, 
		t->expected_heterozygosity(Hₒ_ploidyhet2[1], t, 50),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	V_A2 = map(sum, 2*(0.1)^2 .*island_p2.het)
	plot!(V_A2, grid=false, color=:red, label="\$V_A(t)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 7890f13a-802a-405f-b1b1-ca9857d85be8
plot(p21,p22,p23)

# ╔═╡ 4b791098-5193-46a9-ade3-dd76a630bfaf
t = IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],10), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], β=.25, θ=15.5)

# ╔═╡ 7b2d4c6e-2cd1-4fbd-ae21-c8d87c226253
t.agents

# ╔═╡ 8b559e3b-24f3-4509-abb9-067f70cab4db
fitnesses = exp.(directional_selection(t))


# ╔═╡ 28d218eb-fa98-4ba2-9264-5b1ed9c5debd
truncation(fitnesses, 3)

# ╔═╡ f54e7836-b399-438c-b6af-34f1b7e1296b
island_s2i = IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15., β=0.25)

# ╔═╡ 0f98053e-ec11-4809-9f2e-3e24ca8d26ec
s2i = map(x->(evolving_islanddeme(IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0., θ=15., β=0.25),0.85,20)).pop[end], 1:1000) 

# ╔═╡ 4340b936-744d-47d5-8219-c717a9aa0990
begin
	histogram(s2i, bins = 25, fillalpha = 0.4, title="Estab of diploid migrant")
	xlabel!("\$popsize\$")
	ylabel!("\$n\$")
end

# ╔═╡ a38fbcfc-139c-4642-8226-63c93017309f
length(s2i[s2i .> 0])/1000

# ╔═╡ 0d1dac78-e4b4-4e0f-b9b9-7b1460aca77a
function evolving_islandwbreak(d::IslandDeme, T, M, L, ngen; 
	heterozygosities_p=heterozygosities_p, fit=directional_selection, trait_mean = trait_mean, allelefreqs_p = 
	allelefreqs_p, pf = ploidy_freq, fta = f_trait_agents)
	het = [heterozygosities_p(d)]
	pop = [length(d)]
	tm = [trait_mean(d)]
	af = [allelefreqs_p(d)]
	p2 = [ploidy_freq(d)[2]]
	p3 = [ploidy_freq(d)[3]]
	p4 = [ploidy_freq(d)[4]]
	fta = [f_trait_agents(d)]
	c = 0
	
	for t=1:ngen
		if popsize(d) < L
		migrants = rand(Poisson(M))
		for m in 1:migrants
			migrant = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1)[1]
			push!(d.agents,migrant)
		end
		d = mating_PnB(d,T)
		#d = mutate(d) 
		push!(het, heterozygosities_p(d))
		push!(pop, length(d))
		push!(tm, trait_mean(d))
		push!(af, allelefreqs_p(d))
		push!(p2, ploidy_freq(d)[2])
		push!(p3, ploidy_freq(d)[3])
		push!(p4, ploidy_freq(d)[4])
		push!(fta, f_trait_agents(d))
		c += 1
		else break
		end
	end
	(pop=pop, deme=d, p2=p2, p3=p3, p4=p4, ngen=c, het=het,tm=tm, af=af, fta=fta) 
end

# ╔═╡ c429d513-ec3f-438d-8f01-827ff38befbe
function grid_search(t)
	ploidy = []
	param = []
	pop_size = []
	p2 = []
	p4 = []
	for u in range(0, stop=0.5, length=t)
		for rep in 1:10
		UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.]
		d_p = IslandDeme(agents = randagent_p(0.5, 0.1, 250, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UG, μ=0., θ=12.5, β=0.25)
		sim_ploidyvar = evolving_islandwbreak(d_p, 0.5, 1., 100, 10000)
		if sim_ploidyvar.p2[end] >= sim_ploidyvar.p4[end]
			push!(ploidy,2)
		else
			push!(ploidy,4)
		end
		push!(param, u)
		pop = (sim_ploidyvar.p2[end] + sim_ploidyvar.p4[end])
		push!(pop_size, pop)
		push!(p2, sim_ploidyvar.p2[end])
		push!(p4, sim_ploidyvar.p4[end])
		end
	end
	ploidy, param, pop_size, p2, p4
end

# ╔═╡ 4318b983-7ea8-405e-b62f-c0adf60dd95c
stats = grid_search(100)

# ╔═╡ 196936d0-70ad-4cdf-a248-cf777b2459e1
dp_2 = [(stats[2][i],stats[1][i],stats[3][i]) for i in 1:1000]

# ╔═╡ 3b4cdcc5-fbd8-476b-9a8f-89ca3516de47
begin
function prob(b)
		c = 0
		for x in b
			if x == 4
				c += 1
			end
		end
		c/length(b)
	end

function stabprob(a)
	i = 1
	j = 10
	p = []
	while j <= length(a)
		push!(p,prob(a[i:j]))
		i += 10
		j += 10
	end
	p
	end	
end	

# ╔═╡ 4c4e5667-8893-4ce1-9353-9ff9ede6f818
begin
tick1 = stabprob(stats[1])
plot([0.005:0.005:0.5...],tick1,label=false)
vline!([0.17],label="u=0.17",linewidth=2,style=:dash)
hline!([0.50],label=false,linewidth=2,style=:dash)
end

# ╔═╡ 1bf4ed30-a7cf-42d0-92ec-3aa9369499d9
begin
p3ug = scatter(stats[2],stats[1], label=false, title="UG for 2n and 4n")
vline!([0.17], label="u=0.17",linewidth=5)
xlabel!("u (unreduced gamete formation 2n)")
ylabel!("Ploidy")
end

# ╔═╡ c3f7b7c1-a698-434f-a217-d61338149f05
function bin(data)
	up2 = []
	up4 = []
	for i in data
		if i[3] != 0
			if i[2] == 2
				push!(up2,i[1])
			else
				push!(up4,i[1])
			end
		end
	end
	up2,up4
end

# ╔═╡ 593b3221-8e8c-48a6-9eac-e6e4941a4f1a
binz = bin(dp_2)

# ╔═╡ 4dc218bf-d766-4611-9ddb-f89a3c7b2ad2
begin
	p1ug = histogram(binz[1], bins=50, label=false, title="Diploid > Tetraploid")
	vline!([0.17], label="u=0.17",linewidth=5)
	xlabel!("u (unreduced gamete formation 2n)")
	ylabel!("Count")
end

# ╔═╡ da6ca77e-59a7-4b6d-be7f-f745ddb92efd
begin
	p2ug = histogram(binz[2], bins=50, label=false, title="Tetraploid > Diploid")
	vline!([0.17], label="u=0.17",linewidth=5)
	xlabel!("u (unreduced gamete formation 2n)")
	ylabel!("Count")
end

# ╔═╡ 00239b4e-8087-461d-a78c-54653f645ed0
begin
p7 = scatter(stats[2], stats[3], grid=false, color=:black, label="Pop size after t generations")
scatter!(stats[2], stats[4], grid=false, color=:green, label="Diploids")
scatter!(stats[2], stats[5], grid=false, color=:red, label="Tetraploids")
vline!([0.17], label="u=0.17",linewidth=5)
xlabel!("\$u\$")
ylabel!("Number of individuals")
end

# ╔═╡ Cell order:
# ╟─c016ac10-b047-11eb-37a1-f198dc39687e
# ╟─985ddc8f-d065-449b-9eeb-41003b7449a8
# ╟─a1df4b9d-2ab3-4d20-9399-778330ed200c
# ╠═748df908-da3d-4025-8a0c-369678f82bae
# ╠═7b60ca0b-9ae9-4b41-bd6f-744372135a7a
# ╠═d9f6864f-fdb6-43b8-b8bb-3ea98af682c7
# ╠═11a96686-826f-4875-aec4-2f7d619bc133
# ╠═15902758-b8bd-4b9d-b0d9-46e92007972c
# ╠═d98d3d38-2352-4800-89f4-b077b819d040
# ╠═2342c60c-7e5a-47a3-b23d-a8b8fca19d7c
# ╠═ded1736f-e18e-4174-9a6b-7df33b804c9e
# ╠═2eaa35ba-f42e-4417-82d9-7dc006070d8c
# ╠═a7ce43d2-9ccd-4813-bcd5-b405eb5c0008
# ╠═83cc65c3-3a5f-47fc-8055-7e319000f7da
# ╠═7890f13a-802a-405f-b1b1-ca9857d85be8
# ╠═4b791098-5193-46a9-ade3-dd76a630bfaf
# ╠═7b2d4c6e-2cd1-4fbd-ae21-c8d87c226253
# ╠═8b559e3b-24f3-4509-abb9-067f70cab4db
# ╠═28d218eb-fa98-4ba2-9264-5b1ed9c5debd
# ╠═f54e7836-b399-438c-b6af-34f1b7e1296b
# ╠═0f98053e-ec11-4809-9f2e-3e24ca8d26ec
# ╠═4340b936-744d-47d5-8219-c717a9aa0990
# ╠═a38fbcfc-139c-4642-8226-63c93017309f
# ╠═0d1dac78-e4b4-4e0f-b9b9-7b1460aca77a
# ╠═c429d513-ec3f-438d-8f01-827ff38befbe
# ╠═4318b983-7ea8-405e-b62f-c0adf60dd95c
# ╠═196936d0-70ad-4cdf-a248-cf777b2459e1
# ╟─3b4cdcc5-fbd8-476b-9a8f-89ca3516de47
# ╠═4c4e5667-8893-4ce1-9353-9ff9ede6f818
# ╠═1bf4ed30-a7cf-42d0-92ec-3aa9369499d9
# ╟─c3f7b7c1-a698-434f-a217-d61338149f05
# ╠═593b3221-8e8c-48a6-9eac-e6e4941a4f1a
# ╠═4dc218bf-d766-4611-9ddb-f89a3c7b2ad2
# ╠═da6ca77e-59a7-4b6d-be7f-f745ddb92efd
# ╠═00239b4e-8087-461d-a78c-54653f645ed0
