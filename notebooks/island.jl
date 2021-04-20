### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 1f791df0-8b5e-11eb-3bc5-e7302263af09
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes

# ╔═╡ 325d63e0-8b5e-11eb-0d55-f149adc308e1
using PolyStab

# ╔═╡ 423cffe0-8b5f-11eb-1df4-ff7eed10289a
using PolyStab: Agent, randagent_p, MixedPloidyDeme, IslandDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, random_mating, evolving_neutraldeme, evolving_islanddeme

# ╔═╡ 10d563a1-fb0f-4b6c-ae10-b16f6bef2f86
begin
	using PolyStab:AbstractDeme
	popsize(d::AbstractDeme) = length(d.agents)
end

# ╔═╡ 679ebe33-4d72-4cd8-81d4-bc74fc7bcf60
using PolyStab: mating_PnB, mutate, ploidy_freq

# ╔═╡ a7bc4e02-8b5c-11eb-0e36-ef47de093003
md""" ### Island model
"""

# ╔═╡ ba38834e-8b5c-11eb-246c-99429ab2ae27
md""" ###### cfr. Establishment in a new habitat by polygenic adaptation (Barton and Etheridge, 2017)"""

# ╔═╡ baffcaf2-8b5c-11eb-1479-e110e25be9d5
md""" Representation: large abstract mainland population (assumes: infinite population size in HWLE), migration from mainland to (different smaller) island(s) with finite population size (can be respresented as demes).

Some quick thoughts:
- mainland -> island: establishment of polyploid in a diploid island population
- migration between islands (demes) or only from mainland to island?
- starting population mainland (Gv), island
- colonization event (selfing versus non-selfing, i.e. with nonselfing you need at least two individuals to succesfully colonize an island). 
- stochasticity < dispersal (phenotype - island optimum mismatch), i.e. do all islands have the same optimum as mainland (agents will already be adapted) or is there variation in island optima (different islands with different parameters versus all islands same parameters and different simulation runs).
"""

# ╔═╡ 164cbf80-8b5d-11eb-129c-2fd4863912df
md""" ###### Some simulations"""

# ╔═╡ 3712a3fe-8b5d-11eb-2560-5940bd21ea0c
md""" First we can look a the case with a single migrant migrating to a single island. Simulations start with an idealized infinite mainland population in HWLE. There is migration from this population to the island with migration rate M. Selfing is allowed so a single migrant can establish a new population on an island. The source is assumed to be poorly adapted to conditions on the island with selection gradient β."""

# ╔═╡ d97a6a70-8b5d-11eb-0ed6-056d834ce9bb
md""" The island can be modelled as a single deme:"""

# ╔═╡ 2625fd72-8b5f-11eb-1e69-adf12cbd84c9
island = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0.)

# ╔═╡ c4067530-8fcc-11eb-322e-bdeba99f75a4
islanddeme = IslandDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], β=.2)

# ╔═╡ 8272fcde-8b5f-11eb-1c13-b591229f902f
migrant = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1)[1]

# ╔═╡ 93754ff5-0428-4188-bf46-49137ee1fc4f
trait(migrant)

# ╔═╡ faa45f90-8b6b-11eb-2826-776161a15a58
migrant.loci

# ╔═╡ ce694600-8bd6-11eb-0d50-c7bd17294608
allelefreqs_p(migrant)

# ╔═╡ 05a7f6f0-8bd9-11eb-0892-f3cde611f9ee
mean(allelefreqs_p(migrant))

# ╔═╡ 143183d0-8b61-11eb-1986-310271688189
push!(island.agents,migrant)

# ╔═╡ fc0f65e0-8fcc-11eb-21ab-dd9a09e00355
push!(islanddeme.agents,migrant)

# ╔═╡ 4d1ba2c0-8bcf-11eb-2b00-adaff59d3e9e
length(island.agents)

# ╔═╡ 7b190320-8bcf-11eb-15b6-2b3a871bc3bc
allelefreqs_p(island)

# ╔═╡ 93164980-8bcd-11eb-02a1-496109bc20c3
2 .*heterozygosities_p(island)

# ╔═╡ 91834c70-8bd8-11eb-0ef6-e79beaa1849b
mean(2 .*heterozygosities_p(island))

# ╔═╡ 52c76560-8b61-11eb-1edd-e71d92ca9a1e
island

# ╔═╡ 08db9cd0-8fcd-11eb-1e0a-fd10a73b4dd0
islanddeme

# ╔═╡ 58605950-8b61-11eb-16d9-0ff9be46d489
island_p = evolving_selectiondeme(island,20)

# ╔═╡ 145aba50-8fcd-11eb-081d-4ffe44449e3f
island_p2 = evolving_islanddeme(islanddeme,20)

# ╔═╡ 6df662b0-8b6a-11eb-0e74-611b7747413e
begin
	pf2_p1 = island_p.p2
	pf3_p1 = island_p.p3
	pf4_p1 = island_p.p4
	p1 = plot(pf2_p1, grid=false, color=:blue, label="diploids", legend=:bottomright)
	#plot!(pf3_p1, grid=false, color=:green, label="triploids")
	plot!(pf4_p1, grid=false, color=:red, label="tetraploids")
	hline!([island.K],label ="K",colour = "black",linestyle=:dash)
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ dc3ad710-8b6a-11eb-3e55-7ba4760f8132
begin
	traitmean_ploidy = map(mean, island_p.tm)
	p2 = plot(traitmean_ploidy, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright) #title="Stabilizing selection"
	for (i,t) in enumerate(island_p.fta)
	scatter!([i for x in 1:length(t)],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([island.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ 84585d6e-8c09-11eb-3c5c-af38eb556b15
var_add(a::Agent,α) = ploidy(a)*α^2*sum(heterozygosities_p(a))

# ╔═╡ 7473dcf0-8fcd-11eb-16e8-a5371cf5105f
begin
	pf2_p2 = island_p2.p2
	pf3_p2 = island_p2.p3
	pf4_p2 = island_p2.p4
	p21 = plot(pf2_p2, grid=false, color=:blue, label="diploids", legend=:bottomright)
	#plot!(pf3_p1, grid=false, color=:green, label="triploids")
	plot!(pf4_p2, grid=false, color=:red, label="tetraploids")
	hline!([island.K],label ="K",colour = "black",linestyle=:dash)
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ 1bbe42d0-8fd2-11eb-0d53-8be7e42d7722
begin
	traitmean_ploidy2 = map(mean, island_p2.tm)
	p22 = plot(traitmean_ploidy2, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright) #title="Stabilizing selection"
	for (i,t) in enumerate(island_p2.fta)
	scatter!([i for x in 1:length(t)],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([island.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ 1cd88c70-8fd2-11eb-30bc-a3b3c916b97d
begin
	Hₒ_ploidyhet2 = map(mean, 2 .*island_p2.het)
	expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	p23 = plot(Hₒ_ploidyhet2, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH")
	plot!(1:island_p2.ngen, 
		t->expected_heterozygosity(Hₒ_ploidyhet2[1], t, 50),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	V_A2 = map(sum, 2*(0.1)^2 .*island_p2.het)
	plot!(V_A2, grid=false, color=:red, label="\$V_A(t)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ ee422a80-8b6a-11eb-13ce-a7f8ab3de0ef
begin
	Hₒ_ploidyhet = map(mean, 2 .*island_p.het)
	p3 = plot(Hₒ_ploidyhet, grid=false, color=:black, label="\$H_o(t)\$", title = "LOH")
	plot!(1:island_p.ngen, 
		t->expected_heterozygosity(Hₒ_ploidyhet[1], t, 50),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	V_A = map(sum, 2*(0.1)^2 .*island_p.het)
	plot!(V_A, grid=false, color=:red, label="\$V_A(t)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 5a0f4c50-8c26-11eb-0f10-dfd79a2ff7c5
begin
	plot(p1,p2,p3)
	#savefig("single_migrant.pdf")
end

# ╔═╡ 4be009d0-8fd2-11eb-2338-15f33ab00a20
plot(p21,p22,p23)

# ╔═╡ 4ce3fdbb-82f2-46da-a226-485306b0d91b


# ╔═╡ 602cebcc-46bf-48fe-a3d3-163318796b98
md""" ### Simulations"""

# ╔═╡ 1182e78a-f422-4e14-b6f4-e6f2c05c3da4
island_s2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0.)

# ╔═╡ c3814ccf-976e-4c9c-99dc-ddab6b4945b8
migrant_s2 = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],1)[1]

# ╔═╡ 9db991e8-49af-44be-abdb-646d5409d879
trait(migrant_s2)

# ╔═╡ 6b38a600-768f-4abc-9e3f-63a010b15461
push!(island_s2.agents,migrant_s2)

# ╔═╡ 8841a695-1110-487e-bd8a-720ef80a64c3
s2 = map(x->(evolving_selectiondeme(island_s2,20)).pop[end], 1:1000) 

# ╔═╡ b3f9e34a-2a82-4223-ac2a-6c9e8442b0ab
begin
	histogram(s2, bins = 25, fillalpha = 0.4, title="Estab of diploid migrant")
	xlabel!("\$popsize\$")
	ylabel!("\$n\$")
end

# ╔═╡ 1265d278-030c-48f9-8f2b-b73375b64823
length(s2[s2 .> 0])/1000

# ╔═╡ 06a553ee-662d-4bbf-8a5f-365912b14175
island_s4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 0., 0., 1.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], μ=0.)

# ╔═╡ d7dae018-84ac-493d-96a3-d400073b29a2
migrant_s4 = randagent_p(0.5, 0.5, 50, [0., 0., 0., 1.],1)[1]

# ╔═╡ 2d4c5f74-07a5-4082-b75b-93a772376173
trait(migrant_s4)

# ╔═╡ 65cbc68c-7176-4bc1-8a2a-20f245311918
push!(island_s4.agents,migrant_s4)

# ╔═╡ 0d4614ce-3d3b-4b57-ba53-9b42a11c711b
s4 = map(x->(evolving_selectiondeme(island_s4,20)).pop[end], 1:1000) 

# ╔═╡ 091ec75c-9d54-4b3d-bc45-dbb748ca447a
begin
	histogram(s4, bins = 25, fillalpha = 0.4, title="Estab of tetraploid migrant")
	xlabel!("\$popsize\$")
	ylabel!("\$n\$")
end

# ╔═╡ 558da9c6-e3f1-481e-86c2-8b825ae1f75f
length(s4[s4 .> 0])/1000

# ╔═╡ f29e0d81-0e87-4d39-adaa-57dc080497c9
md""" oooo """

# ╔═╡ 4abbc1ed-6bdf-4bf0-8cbd-3deee6770543
md""" ##### Should be straightforward to implement habitat in 2D"""

# ╔═╡ 4b626564-9481-475d-a9f5-535ad7edc11f
begin
d0 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.85 0.15 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], K=50)
d10 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],10), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], K=50)
d50 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],50), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.85 0.15 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], K=50)
d500 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],500), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.85 0.15 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], K=50)
end

# ╔═╡ f64591c7-ee98-49e9-b849-25528c367757
h_2D = [d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d50 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0; d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0]

# ╔═╡ 5605ac82-c3f1-4869-b82e-70c8e6c5ff73
begin
	data = rand(10,10)
	heatmap(data)
end

# ╔═╡ 493ea7a3-c639-4519-8288-df558cf6debf
@with_kw struct Habitat2D{D,T}
    demes::Matrix{D}
    σ ::T = sqrt(1/2) 
    b ::T = 0.1 
	θ ::T = 12.5 
    Dm::T = 250. 
end

# ╔═╡ 09972d7b-61fa-46d3-86c3-abe6084334c8
begin
hab2D = Habitat2D(demes=h_2D)
popsize.(hab2D.demes)
heatmap(popsize.(h_2D), c=cgrad([:white,:red,:yellow]))
end

# ╔═╡ 5afa6717-4b0c-481b-9d8c-66f37ab977d7
#emptycopy(h::Habitat2D) = Habitat2D(emptycopy.(h.demes), h.σ, h.b, h.θ, h.Dm)

# ╔═╡ 601d4159-a431-48e7-925a-d9fa2e0f354f
emptycopy(d::MixedPloidyDeme{A,T}) where A where T = 
    MixedPloidyDeme(A[], d.K, d.θ, d.rm, d.Vs, d.α, d.μ, d.OV, d.UG)

# ╔═╡ ec7d7863-f16b-4667-87fd-f32383b7ced0
function emptycopy(h::Habitat2D)
	new_h = similar(h.demes)
	for i in 1:size(h.demes)[1]
		for j in 1:size(h.demes)[2]
			new_h[i,j] = emptycopy(h.demes[i,j])
		end
	end
	Habitat2D(demes=new_h, σ=h.σ, b=h.b, θ=h.θ, Dm=h.Dm)
end

# ╔═╡ 43cca8f0-a2f3-4336-84d0-edd96910596a
function random_walk2D(h::Habitat2D, p)
    new_h = emptycopy(h)
	for (i, deme) in pairs(h.demes)
        for agent in deme.agents
            step = (rand(),rand()) < (p,p) ? rand([(-1,0),(1,0),(0,-1),(0,1),(-1,1),(-1,-1),(1,-1),(1,1)]) : (0,0) 
            if step[1] == -1 && i[1] == 1
                step = (0,0)
            elseif step[1] == 1  && i[1] == size(h.demes)[1]
                step = (0,0)
			elseif step[2] == -1  && i[2] == 1
                step = (0,0)
			elseif step[2] == 1  && i[2] == size(h.demes)[2]
                step = (0,0)
            end
            push!(new_h.demes[i[1]+step[1],i[2]+step[2]], agent)
        end
    end
    new_h
end

# ╔═╡ 332e259d-c58f-4c70-b4b6-f6a5d6f8bfb9
function gaussian_dispersal2D(h::Habitat2D,σ)
    new_h = emptycopy(h)
	dist = Normal(0,σ)
	dist_trunc = truncated(dist,-2*σ,2*σ)
	bin_1 = pdf(dist, σ)
	for (i, deme) in pairs(h.demes)
        for agent in deme.agents
            step = -bin_1 < rand(dist_trunc) < bin_1  ?  (0,0) : rand([(-1,0),(1,0),(0,-1),(0,1),(-1,1),(-1,-1),(1,-1),(1,1)])
            if step[1] == -1 && i[1] == 1
                step = (0,0)
            elseif step[1] == 1  && i[1] == size(h.demes)[1]
                step = (0,0)
			elseif step[2] == -1  && i[2] == 1
                step = (0,0)
			elseif step[2] == 1  && i[2] == size(h.demes)[2]
                step = (0,0)
            end
            push!(new_h.demes[i[1]+step[1],i[2]+step[2]], agent)
        end
    end
    new_h
end

# ╔═╡ 2cd2d5d7-1be2-4d21-a0a2-21dc2957e36b
new_h = emptycopy(hab2D)

# ╔═╡ 2ab05bac-1d9e-4e87-bd38-2797bde88305
emptycopy(hab2D)

# ╔═╡ 994f4457-7425-4b43-a12a-c00dd25afe9f
function random_walk2D(h::Habitat2D,σ,ngen)
	data = []
	for n = 1:ngen
		h = random_walk2D(h,σ)
		pop = popsize.(h.demes)
		push!(data, pop)
	end
	(h=h, ngen=ngen, data=data)
end
	

# ╔═╡ b58b39fc-e454-4a98-8952-9e9b736ff360
function gaussian_dispersal2D(h::Habitat2D,σ,ngen)
	data = []
	for n = 1:ngen
		h = gaussian_dispersal2D(h,σ)
		pop = popsize.(h.demes)
		push!(data, pop)
	end
	(h=h, ngen=ngen, data=data)
end
	

# ╔═╡ 66d4e853-aa19-4c35-b8ca-5592e0f211c1
begin
hab2Dwt = gaussian_dispersal2D(hab2D,0.5)
popsize.(hab2Dwt.demes)
heatmap(popsize.(hab2Dwt.demes), c=cgrad([:white,:red,:yellow]))
end

# ╔═╡ 32a50950-9d75-4b23-b27e-8c9f19b3a5dd
sim2Dgaus = gaussian_dispersal2D(hab2D,0.5,200)

# ╔═╡ be9268b2-18fa-43cc-b266-fa3305149865
#sim2Drand = random_walk2D(hab2D,0.2,200)

# ╔═╡ dfbc5a08-ea5d-4e48-9f10-765914ff8bd7
begin
	anim_ranger = @animate for i ∈ 1:250
		heatmap(sim2Drand.data[i], c=cgrad([:white,:grey,:black]))
	end every 1
	gif(anim_ranger, "popsize2Drand.gif", fps = 5)
end

# ╔═╡ 8aea6c5b-b06a-4095-859c-515e2e21ceb6
begin
	anim_rangeg = @animate for i ∈ 1:200
		heatmap(sim2Dgaus.data[i], c=cgrad([:white,:grey,:black]))
	end every 1
	gif(anim_rangeg, "popsize2Dgaus.gif", fps = 4)
end

# ╔═╡ 1fcba851-f04d-4569-9b47-439511b7371f
function evolving_habitat2D(h::Habitat2D{D}, ngen) where D
	data = []
	for n = 1:ngen
		h = random_walk2D(h,0.5)
		#ih = Gaussian_dispersal(h,σ)
		new_h = Matrix{D}(undef, size(h.demes))
		for (i, d) in pairs(h.demes)
			d = mating_PnB(d)
			d = mutate(d)
			new_h[i[1],i[2]] = d
		end
		h = Habitat2D(demes=new_h)
		push!(data, h)
	end
	(h=h, ngen=ngen, data=data)
end

# ╔═╡ 990052c3-83a2-4101-b46c-1a6bccd91a5b
sim_hab = evolving_habitat2D(hab2D, 15)

# ╔═╡ ec7ca26d-7854-4144-b8b8-cbe647c02f1e
begin
get_first(v)=v[1]
get_second(v)=v[2]
get_third(v)=v[3]
get_fourth(v)=v[4]
end

# ╔═╡ 013f1d53-9d90-46d4-ad82-aa22115f7435
begin
	K = d50.K
	σ = hab2D.σ
	b = hab2D.b
	Vs = d50.Vs
	rm = d50.rm
	Dm = hab2D.Dm
	s = 1
	
	pop_sizes = popsize.(sim_hab[1].demes)
	ppf1 = get_second.(ploidy_freq.(sim_hab[1].demes))
	ppf2 = get_third.(ploidy_freq.(sim_hab[1].demes))
	ppf3 = get_fourth.(ploidy_freq.(sim_hab[1].demes))
	
	p1h2d = heatmap(pop_sizes, c=cgrad([:white,:grey,:black]))

	p2h2d = heatmap(ppf1, c=cgrad([:white,:blue]))
	p3h2d = heatmap(ppf2, c=cgrad([:white,:green]))
	p4h2d = heatmap(ppf3, c=cgrad([:white,:red]))
	
	

end

# ╔═╡ 57d5b41e-9c5b-4107-8c08-6545d6d3609c
plot(p1h2d,p2h2d,p3h2d,p4h2d)

# ╔═╡ Cell order:
# ╟─a7bc4e02-8b5c-11eb-0e36-ef47de093003
# ╟─ba38834e-8b5c-11eb-246c-99429ab2ae27
# ╟─baffcaf2-8b5c-11eb-1479-e110e25be9d5
# ╟─164cbf80-8b5d-11eb-129c-2fd4863912df
# ╠═1f791df0-8b5e-11eb-3bc5-e7302263af09
# ╠═325d63e0-8b5e-11eb-0d55-f149adc308e1
# ╠═423cffe0-8b5f-11eb-1df4-ff7eed10289a
# ╟─3712a3fe-8b5d-11eb-2560-5940bd21ea0c
# ╟─d97a6a70-8b5d-11eb-0ed6-056d834ce9bb
# ╠═2625fd72-8b5f-11eb-1e69-adf12cbd84c9
# ╠═c4067530-8fcc-11eb-322e-bdeba99f75a4
# ╠═8272fcde-8b5f-11eb-1c13-b591229f902f
# ╠═93754ff5-0428-4188-bf46-49137ee1fc4f
# ╠═faa45f90-8b6b-11eb-2826-776161a15a58
# ╠═ce694600-8bd6-11eb-0d50-c7bd17294608
# ╠═05a7f6f0-8bd9-11eb-0892-f3cde611f9ee
# ╠═143183d0-8b61-11eb-1986-310271688189
# ╠═fc0f65e0-8fcc-11eb-21ab-dd9a09e00355
# ╠═4d1ba2c0-8bcf-11eb-2b00-adaff59d3e9e
# ╠═7b190320-8bcf-11eb-15b6-2b3a871bc3bc
# ╠═93164980-8bcd-11eb-02a1-496109bc20c3
# ╠═91834c70-8bd8-11eb-0ef6-e79beaa1849b
# ╠═52c76560-8b61-11eb-1edd-e71d92ca9a1e
# ╠═08db9cd0-8fcd-11eb-1e0a-fd10a73b4dd0
# ╠═58605950-8b61-11eb-16d9-0ff9be46d489
# ╠═145aba50-8fcd-11eb-081d-4ffe44449e3f
# ╠═6df662b0-8b6a-11eb-0e74-611b7747413e
# ╠═dc3ad710-8b6a-11eb-3e55-7ba4760f8132
# ╠═ee422a80-8b6a-11eb-13ce-a7f8ab3de0ef
# ╠═84585d6e-8c09-11eb-3c5c-af38eb556b15
# ╠═7473dcf0-8fcd-11eb-16e8-a5371cf5105f
# ╠═1bbe42d0-8fd2-11eb-0d53-8be7e42d7722
# ╠═1cd88c70-8fd2-11eb-30bc-a3b3c916b97d
# ╠═5a0f4c50-8c26-11eb-0f10-dfd79a2ff7c5
# ╠═4be009d0-8fd2-11eb-2338-15f33ab00a20
# ╠═4ce3fdbb-82f2-46da-a226-485306b0d91b
# ╟─602cebcc-46bf-48fe-a3d3-163318796b98
# ╠═1182e78a-f422-4e14-b6f4-e6f2c05c3da4
# ╠═c3814ccf-976e-4c9c-99dc-ddab6b4945b8
# ╠═9db991e8-49af-44be-abdb-646d5409d879
# ╠═6b38a600-768f-4abc-9e3f-63a010b15461
# ╠═8841a695-1110-487e-bd8a-720ef80a64c3
# ╠═b3f9e34a-2a82-4223-ac2a-6c9e8442b0ab
# ╠═1265d278-030c-48f9-8f2b-b73375b64823
# ╠═06a553ee-662d-4bbf-8a5f-365912b14175
# ╠═d7dae018-84ac-493d-96a3-d400073b29a2
# ╠═2d4c5f74-07a5-4082-b75b-93a772376173
# ╠═65cbc68c-7176-4bc1-8a2a-20f245311918
# ╠═0d4614ce-3d3b-4b57-ba53-9b42a11c711b
# ╠═091ec75c-9d54-4b3d-bc45-dbb748ca447a
# ╠═558da9c6-e3f1-481e-86c2-8b825ae1f75f
# ╠═f29e0d81-0e87-4d39-adaa-57dc080497c9
# ╟─4abbc1ed-6bdf-4bf0-8cbd-3deee6770543
# ╠═4b626564-9481-475d-a9f5-535ad7edc11f
# ╠═f64591c7-ee98-49e9-b849-25528c367757
# ╠═10d563a1-fb0f-4b6c-ae10-b16f6bef2f86
# ╠═5605ac82-c3f1-4869-b82e-70c8e6c5ff73
# ╠═493ea7a3-c639-4519-8288-df558cf6debf
# ╠═09972d7b-61fa-46d3-86c3-abe6084334c8
# ╠═5afa6717-4b0c-481b-9d8c-66f37ab977d7
# ╠═601d4159-a431-48e7-925a-d9fa2e0f354f
# ╠═ec7d7863-f16b-4667-87fd-f32383b7ced0
# ╠═43cca8f0-a2f3-4336-84d0-edd96910596a
# ╠═332e259d-c58f-4c70-b4b6-f6a5d6f8bfb9
# ╠═2cd2d5d7-1be2-4d21-a0a2-21dc2957e36b
# ╠═2ab05bac-1d9e-4e87-bd38-2797bde88305
# ╠═66d4e853-aa19-4c35-b8ca-5592e0f211c1
# ╠═994f4457-7425-4b43-a12a-c00dd25afe9f
# ╠═b58b39fc-e454-4a98-8952-9e9b736ff360
# ╠═32a50950-9d75-4b23-b27e-8c9f19b3a5dd
# ╠═be9268b2-18fa-43cc-b266-fa3305149865
# ╠═dfbc5a08-ea5d-4e48-9f10-765914ff8bd7
# ╠═8aea6c5b-b06a-4095-859c-515e2e21ceb6
# ╠═679ebe33-4d72-4cd8-81d4-bc74fc7bcf60
# ╠═1fcba851-f04d-4569-9b47-439511b7371f
# ╠═990052c3-83a2-4101-b46c-1a6bccd91a5b
# ╠═013f1d53-9d90-46d4-ad82-aa22115f7435
# ╠═ec7ca26d-7854-4144-b8b8-cbe647c02f1e
# ╠═57d5b41e-9c5b-4107-8c08-6545d6d3609c
