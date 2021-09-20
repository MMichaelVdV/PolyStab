### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 62e605b0-a299-11eb-2e86-6fa01314f1ec
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes

# ╔═╡ 8147adb8-241e-4eaa-9e7d-333281d91021
using PolyStab: Agent, randagent_p, MixedPloidyDeme, IslandDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, random_mating, evolving_neutraldeme, evolving_islanddeme

# ╔═╡ 1b1b41f7-6549-4e57-8a40-ea06d3f0d957
using PolyStab: mating_PnB, mutate, ploidy_freq

# ╔═╡ 6bf63cc4-ed68-4e4b-a2ad-ed19e3d8158b
begin
	using PolyStab:AbstractDeme
	popsize(d::AbstractDeme) = length(d.agents)
end

# ╔═╡ 0513c900-a191-4f1a-80c4-c3c4045a0d6c
using Measures

# ╔═╡ fe349ec8-3a44-44d0-b131-1fbe8224b2ae
using PolyStab: trait_mean

# ╔═╡ a449125a-0d09-41c2-ab28-1afd40065747
begin
d0 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ;  0. 1. 0. 0. ], K=50)
d0p = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ;  0. 1. 0. 0. ], K=50, θ=13.5)
d10 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],10), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ;  0. 1. 0. 0. ], K=50)
d50 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],50), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ;  0. 1. 0. 0. ], K=50)
d500 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],500), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0. ], K=50)
end

# ╔═╡ 45aa806a-c591-4c27-a5e6-c8cb412561f7
h_2D = [d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p; d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p; d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d50 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0 d0p d0p d0p; d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p; d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p; d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p d0p]

# ╔═╡ a7433596-a76b-405c-a8cb-7ee15ee3ab55
begin
	da = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0.5 0.5 0. 0. ;  0. 1. 0. 0. ], K=50)
	db = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0.5 0.5 0. 0. ;  0. 1. 0. 0. ], K=50, θ=13.)
	dc = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0.5 0.5 0. 0. ;  0. 1. 0. 0. ], K=50, θ=14.)
	dd = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0.5 0.5 0. 0. ;  0. 1. 0. 0. ], K=50, θ=15.)
	de = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0.5 0.5 0. 0. ;  0. 1. 0. 0. ], K=50, θ=16.)
	df = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0.5 0.5 0. 0. ;  0. 1. 0. 0. ], K=50, θ=17.)
end

# ╔═╡ 330a734f-7e8f-415f-9997-ef24b3d4b654
h_2D_rainbow = [da da da da da da da da d50 d50 da da da da da da da da; da da da da da da da da da da da da da da da da da da; db db db db db db db db db db db db db db db db db db; db db db db db db db db db db db db db db db db db db; dc dc dc dc dc dc dc dc dc dc dc dc dc dc dc dc dc dc; dc dc dc dc dc dc dc dc dc dc dc dc dc dc dc dc dc dc; dd dd dd dd dd dd dd dd dd dd dd dd dd dd dd dd dd dd; dd dd dd dd dd dd dd dd dd dd dd dd dd dd dd dd dd dd; de de de de de de de de de de de de de de de de de de; de de de de de de de de de de de de de de de de de de; df df df df df df df df df df df df df df df df df df; df df df df df df df df df df df df df df df df df df]

# ╔═╡ e2249099-6895-4b9f-9b51-79d84d9a4238
pheno(deme)=deme.θ

# ╔═╡ c8b4e37d-f8e4-4ddf-b2eb-a6cc8c2ec71e
begin
heatmap(pheno.(h_2D_rainbow), c=cgrad([:white,:grey,:black]), title="Environmental gradient",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=14, legendfontsize=14,xticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],yticks=[1,2,3,4,5,6,7,8,9,10,11,12])
savefig("envgrad2d")
end

# ╔═╡ f5e391cc-92a2-4918-ab33-61e1b70a0e2b
begin
	data = rand(10,10)
	heatmap(data)
end

# ╔═╡ 1ec006fd-1938-4614-8c23-e4280a245117
@with_kw struct Habitat2D{D,T}
    demes::Matrix{D}
    σ ::T = sqrt(1/2) 
    b ::T = 0.1 
	θ ::T = 12.5 
    Dm::T = 250. 
end

# ╔═╡ b1f0127f-1d83-45bb-90aa-5e1934d6f634
begin
hab2D = Habitat2D(demes=h_2D_rainbow, σ=.1)
popsize.(hab2D.demes)
heatmap(popsize.(h_2D), c=cgrad([:white,:red,:yellow]))
end

# ╔═╡ 2c119302-81be-4d34-a8ea-4d9156acb0e9
emptycopy(d::MixedPloidyDeme{A,T}) where A where T = 
    MixedPloidyDeme(A[], d.K, d.θ, d.rm, d.Vs, d.α, d.μ, d.OV, d.UG)

# ╔═╡ 1d0c83e5-01af-4180-8684-f2b7a99ae2e4
function emptycopy(h::Habitat2D)
	new_h = similar(h.demes)
	for i in 1:size(h.demes)[1]
		for j in 1:size(h.demes)[2]
			new_h[i,j] = emptycopy(h.demes[i,j])
		end
	end
	Habitat2D(demes=new_h, σ=h.σ, b=h.b, θ=h.θ, Dm=h.Dm)
end

# ╔═╡ fac3cf77-1e80-4269-989f-856ff986d807
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

# ╔═╡ cf2cbacc-0ce1-4444-825b-49ca34e4c008
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

# ╔═╡ ab54503f-44c1-43fa-a012-f7a22ab4d531
new_h = emptycopy(hab2D)

# ╔═╡ 0aa4cabc-55f7-4ac1-a09e-212162296de6
emptycopy(hab2D)

# ╔═╡ 4b5d6484-d8b2-447d-9713-52a127fb0a47
function random_walk2D(h::Habitat2D,σ,ngen)
	data = []
	for n = 1:ngen
		h = random_walk2D(h,σ)
		pop = popsize.(h.demes)
		push!(data, pop)
	end
	(h=h, ngen=ngen, data=data)
end

# ╔═╡ 768e7c97-bcbc-4f18-a2cd-7867da681950
function gaussian_dispersal2D(h::Habitat2D,σ,ngen)
	data = []
	for n = 1:ngen
		h = gaussian_dispersal2D(h,σ)
		pop = popsize.(h.demes)
		push!(data, pop)
	end
	(h=h, ngen=ngen, data=data)
end

# ╔═╡ c9bedbbd-5e2c-40fa-b010-15d50d2752d9
begin
hab2Dwt = gaussian_dispersal2D(hab2D,0.5)
popsize.(hab2Dwt.demes)
heatmap(popsize.(hab2Dwt.demes), c=cgrad([:white,:red,:yellow]))
end

# ╔═╡ bacec843-9f78-4360-bd65-c181c3aa88ef
sim2Dgaus = gaussian_dispersal2D(hab2D,0.1,200)

# ╔═╡ 306daeac-8648-4ac4-a2d5-ccae385f517b
begin
	anim_rangeg = @animate for i ∈ 1:200
		heatmap(sim2Dgaus.data[i], c=cgrad([:white,:grey,:black]))
	end every 1
	gif(anim_rangeg, "popsize2Dgaus.gif", fps = 4)
end

# ╔═╡ 2113c2d4-e16c-488c-96c8-6d7ac3a458ce
function evolving_habitat2D(h::Habitat2D{D}, ngen) where D
	data = []
	for n = 1:ngen
		h = random_walk2D(h,0.15)
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

# ╔═╡ 801ae3ca-3e92-4f00-9c15-cb277aa01228
sim_hab = evolving_habitat2D(hab2D, 250)

# ╔═╡ 595d5523-e9f3-4652-9ea9-1b1e9883a13d
begin
get_first(v)=v[1]
get_second(v)=v[2]
get_third(v)=v[3]
get_fourth(v)=v[4]
end

# ╔═╡ 1352e25f-89b2-44d5-a41b-f9e615ca2887
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
	ppf4 = argmax.(ploidy_freq.(sim_hab[1].demes))
	
	ph2d = heatmap(pop_sizes, clim=(0.,65.), c=cgrad([:white,:grey,:black]), title="Population size",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=14, legendfontsize=14,xticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],yticks=[1,2,3,4,5,6,7,8,9,10,11,12])
	
	p1h2d = heatmap(ppf1,clim=(0.,65.), c=cgrad([:white,:blue]),title="Diploid",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=14, legendfontsize=14,xticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],yticks=[1,2,3,4,5,6,7,8,9,10,11,12])
	
	p2h2d = heatmap(ppf2,clim=(0.,65.), c=cgrad([:white,:purple]),title="Triploid",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=14, legendfontsize=14,xticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],yticks=[1,2,3,4,5,6,7,8,9,10,11,12])
	
	p3h2d = heatmap(ppf3,clim=(0.,65.), c=cgrad([:white,:red]),title="Tetraploid",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=14, legendfontsize=14,xticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],yticks=[1,2,3,4,5,6,7,8,9,10,11,12])
	
	p4h2d = heatmap(ppf4, c=cgrad([:white,:blue,:purple,:red], categorical=true), clim=(0.,4.),title="Most abundant ploidy",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=14, legendfontsize=14,xticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],yticks=[1,2,3,4,5,6,7,8,9,10,11,12])
	
end

# ╔═╡ dc70d416-7ac5-4c06-acce-9926ce9fe399
savefig(ph2d,"start2dsim")

# ╔═╡ 66d2cf4c-f6b1-462c-b124-f3eb1dc9adb8
begin
	anim_2D = @animate for i ∈ 1:100
		heatmap(argmax.(ploidy_freq.(sim_hab.data[i].demes)), c=cgrad([:white,:blue,:green,:red]), clim= (-0.1,4.))
		end every 1
	gif(anim_2D, "polystab2Dgaus.gif", fps = 4)
end

# ╔═╡ d7fd219f-7c3d-49ae-8143-e7144d3db253
sim_hab.data

# ╔═╡ 15b51f89-7488-4c02-941e-068bc53ebfe6
heatmap(sim_hab.data)

# ╔═╡ a1d9ce5f-a5fc-4fbd-8841-9f2763be7e06
sim_hab.data

# ╔═╡ 221bf2a9-afb4-4d98-96bc-248960c1f99b
begin
	function f_het_demes(h::Habitat2D{D}) where D
    het_demes = Matrix{Float64}(undef, size(h.demes))
	for (i, d) in pairs(h.demes)
		if length(d) != 0
			het_demes[i[1],i[2]] = mean(heterozygosities_p(d))
		else
			het_demes[i[1],i[2]] = 0.
		end
	end
	het_demes
	end
	
anim_range_Vg = @animate for i ∈ 1:100
	sim_habA = sim_hab.data[i]

	het_demes = f_het_demes(sim_habA)
		
	p1 = heatmap(het_demes, grid=false, c=cgrad([:white,:black]), label="Vg_mean deme")
	
	xlabel!("Space")
	ylabel!("Genetic variance")
	end every 1
	gif(anim_range_Vg, "geneticvariance2D.gif", fps = 3)
end

# ╔═╡ 035059d7-8d84-4077-87c4-a05d7ab3e1cc
p_het = heatmap(f_het_demes(sim_hab.h), grid=false, c=cgrad([:white,:black]), label="Vg_mean deme", title="Genetic variance", xticks=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],yticks=[1,2,3,4,5,6,7,8,9,10,11,12],xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=14, legendfontsize=14)

# ╔═╡ e49cf0fe-7143-434a-bae8-e1103f2a3701
begin
	plot(ph2d,p4h2d,p_het,p1h2d,p2h2d,p3h2d, size=(1600,900), margins=1mm)
	savefig("hab2d")
end

# ╔═╡ a76cce22-c6d3-4863-a674-8b1daa95444a
begin
	function f_trait_demes(h::Habitat2D{D}) where D
    trait_demes = Matrix{Float64}(undef, size(h.demes))
	for (i, d) in pairs(h.demes)
		if length(d) != 0
			trait_demes[i[1],i[2]] = trait_mean(d)
		else
			trait_demes[i[1],i[2]] = 0.
		end
	end
	trait_demes
	end
	
anim_range_tm = @animate for i ∈ 1:150
	sim_habA = sim_hab.data[i]

	trait_demes = f_trait_demes(sim_habA)
		
	p1 = heatmap(trait_demes, grid=false, c=cgrad([:red,:orange,:yellow,:green,:blue,:purple]), label="trait_mean deme")
	
	#xlabel!("Space")
	#ylabel!("Trait mean")
	end every 1
	gif(anim_range_tm, "traitmean2D.gif", fps = 5)
end

# ╔═╡ Cell order:
# ╠═62e605b0-a299-11eb-2e86-6fa01314f1ec
# ╠═8147adb8-241e-4eaa-9e7d-333281d91021
# ╠═1b1b41f7-6549-4e57-8a40-ea06d3f0d957
# ╠═a449125a-0d09-41c2-ab28-1afd40065747
# ╠═45aa806a-c591-4c27-a5e6-c8cb412561f7
# ╠═a7433596-a76b-405c-a8cb-7ee15ee3ab55
# ╠═330a734f-7e8f-415f-9997-ef24b3d4b654
# ╠═e2249099-6895-4b9f-9b51-79d84d9a4238
# ╠═c8b4e37d-f8e4-4ddf-b2eb-a6cc8c2ec71e
# ╠═6bf63cc4-ed68-4e4b-a2ad-ed19e3d8158b
# ╠═f5e391cc-92a2-4918-ab33-61e1b70a0e2b
# ╠═1ec006fd-1938-4614-8c23-e4280a245117
# ╠═b1f0127f-1d83-45bb-90aa-5e1934d6f634
# ╠═2c119302-81be-4d34-a8ea-4d9156acb0e9
# ╠═1d0c83e5-01af-4180-8684-f2b7a99ae2e4
# ╠═fac3cf77-1e80-4269-989f-856ff986d807
# ╠═cf2cbacc-0ce1-4444-825b-49ca34e4c008
# ╠═ab54503f-44c1-43fa-a012-f7a22ab4d531
# ╠═0aa4cabc-55f7-4ac1-a09e-212162296de6
# ╠═c9bedbbd-5e2c-40fa-b010-15d50d2752d9
# ╠═4b5d6484-d8b2-447d-9713-52a127fb0a47
# ╠═768e7c97-bcbc-4f18-a2cd-7867da681950
# ╠═bacec843-9f78-4360-bd65-c181c3aa88ef
# ╠═306daeac-8648-4ac4-a2d5-ccae385f517b
# ╠═2113c2d4-e16c-488c-96c8-6d7ac3a458ce
# ╠═801ae3ca-3e92-4f00-9c15-cb277aa01228
# ╠═1352e25f-89b2-44d5-a41b-f9e615ca2887
# ╠═dc70d416-7ac5-4c06-acce-9926ce9fe399
# ╠═595d5523-e9f3-4652-9ea9-1b1e9883a13d
# ╠═035059d7-8d84-4077-87c4-a05d7ab3e1cc
# ╠═0513c900-a191-4f1a-80c4-c3c4045a0d6c
# ╠═e49cf0fe-7143-434a-bae8-e1103f2a3701
# ╠═66d2cf4c-f6b1-462c-b124-f3eb1dc9adb8
# ╠═d7fd219f-7c3d-49ae-8143-e7144d3db253
# ╠═15b51f89-7488-4c02-941e-068bc53ebfe6
# ╠═a1d9ce5f-a5fc-4fbd-8841-9f2763be7e06
# ╠═221bf2a9-afb4-4d98-96bc-248960c1f99b
# ╠═fe349ec8-3a44-44d0-b131-1fbe8224b2ae
# ╠═a76cce22-c6d3-4863-a674-8b1daa95444a
