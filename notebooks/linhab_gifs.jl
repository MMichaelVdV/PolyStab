### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 80bd1070-8e89-11eb-2348-954f56cd5b7b
using Statistics, Plots 

# ╔═╡ ddff89b0-8e58-11eb-2fef-5b94ca17af89
using PolyStab: randagent_p, MixedPloidyDeme, Habitat, linear_gradient, initiate_habitat, evolving_habitat, ploidy_freq, heterozygosities_p, trait_mean, trait

# ╔═╡ 381b957c-e4d7-46bd-b164-c26ce9f3bf7c
d = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.], 0), K=25,  OV = [1.0 0.1 0.0 0.0; 0.1 1.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], UG = [0.0 0.0 0.0 0.0; 0.92 0.08 0.0 0.0; 0.5 0.5 0. 0.; 0. 1. 0.0 0.], Vs=1.2)

# ╔═╡ 22341570-8e62-11eb-2dc0-919c8bbaeb55
habi = Habitat(demes=[d])

# ╔═╡ 5bc3759e-8f3f-11eb-3458-3d39fe093824
g_lin = linear_gradient(0.1,12.5,250)

# ╔═╡ 5c000970-8f3f-11eb-0c5d-ef59668300a9
hab = initiate_habitat(d, g_lin, 0.5, 0.5, 50, 20)

# ╔═╡ 675ae674-5132-480d-87c5-2e0ad0c6605b
habs = Habitat(demes=hab.demes, σ=hab.σ, b=0.1, θ=12.5, Dm=250.)

# ╔═╡ 8e623f10-8e85-11eb-098f-7184885620f4
sim_hab = evolving_habitat(habs, 300)

# ╔═╡ 9f274e80-8e85-11eb-1b80-f5511b4a63f0
begin
	K = d.K
	σ = hab.σ
	b = hab.b
	Vs = d.Vs
	rm = d.rm
	Dm = hab.Dm
	s = 1
	pop_sizes = [length(deme) for deme  in sim_hab[1].demes]
	margin = (sqrt(2)*b*σ)/((2*rm*sqrt(Vs))-b*σ) .>= 0.15.*pop_sizes.*σ*sqrt(s)
	ppf1 = [ploidy_freq(deme)[2] for deme  in sim_hab[1].demes]
	ppf2 = [ploidy_freq(deme)[3] for deme  in sim_hab[1].demes]
	ppf3 = [ploidy_freq(deme)[4] for deme  in sim_hab[1].demes]
	pop_sizes = map(Statistics.mean, pop_sizes)
	p1 = plot(pop_sizes, grid=false, color=:black, label=false, legendfontsize=5)
	plot!(pop_sizes, grid=false, color=:black, label=false, legendfontsize=5,linestyle=:dash)
	hline!([K], label = "K")
	hline!([K*(1-(σ*b)*(1/(2*sqrt(Vs)*rm)))], label = "Expected pop size")
	vline!([Dm/2], label = "Starting deme")
	plot!([margin]*10, label = "Deterministic range margin")
	plot!(ppf1, grid=false, color=:blue, label="Diploid")
	plot!(ppf2, grid=false, color=:green, label="Triploid")
	plot!(ppf3, grid=false, color=:red, label="Tetraploid")
	xlabel!("Space")
	ylabel!("Population size N")
	savefig(p1, "popsize1D")
end

# ╔═╡ 82bcd140-8efb-11eb-1d4a-cf573d3bffe0
begin
	function f_het_demes(h::Habitat)
	het_demes = []
	cordsh = []
	for (i, deme) in enumerate(h.demes)
	if length(deme) != 0
		het = d.α^2*sum(heterozygosities_p(deme))
		push!(cordsh,i)
		push!(het_demes,het)
	else
		push!(cordsh,i)
		push!(het_demes,0)
		end
	end
	het_demes, cordsh
	end
	anim_range_Vg = @animate for i ∈ 1:300
		sim_habA = sim_hab.data[i]
		#if i != 1
		#sim_habA = evolving_habitat(sim_habA[1],1,1.06,0.5,10^-6,0.50)
		#end
		pop_sizes = [length(deme) for deme  in sim_habA.demes]
		pop_sizes_p = map(mean, pop_sizes)
		
		#./pop_sizes
		
		ppf1 = [ploidy_freq(deme)[2] for deme  in sim_habA.demes]/100
		ppf2 = [ploidy_freq(deme)[3] for deme  in sim_habA.demes]/100
		ppf3 = [ploidy_freq(deme)[4] for deme  in sim_habA.demes]/100
		het_demes, cordsh = f_het_demes(sim_habA)
		
		p1 = plot(cordsh, het_demes, grid=false, color=:black, label="Vg_mean deme")
		vline!([Dm/2], label = "Starting deme")
		#plot!([margin]*1, label = "Deterministic range margin")
		#plot!(ppf1, grid=false, color=:blue, label="Diploid")
		#plot!(ppf2, grid=false, color=:green, label="Triploid")
		#plot!(ppf3, grid=false, color=:red, label="Tetraploid")
		
		
		xlabel!("Space")
		ylabel!("Genetic variance")
	end every 1
	gif(anim_range_Vg, "geneticvariance.gif", fps = 4)
end

# ╔═╡ e34a6b20-8efc-11eb-059d-397e876bd4bc
begin
	anim_range = @animate for i ∈ 1:300
		sim_habA = sim_hab.data[i]
		#if i != 1
		#sim_habA = evolving_habitat(sim_habA[1],1,1.06,0.5,10^-6,0.50)
		#end
		pop_sizes = [length(deme) for deme  in sim_habA.demes]
		pop_sizes_p = map(mean, pop_sizes)
		
		ppf1 = [ploidy_freq(deme)[2] for deme  in sim_habA.demes]
		ppf2 = [ploidy_freq(deme)[3] for deme  in sim_habA.demes]
		ppf3 = [ploidy_freq(deme)[4] for deme  in sim_habA.demes]
		het_demes, cordsh = f_het_demes(sim_habA)
		
		#plot(pop_sizes_p, grid=false, color=:black, label="Tetraploid")
		#hline!([K], label = "K")
		#hline!([K*(1-(σ*b)*(1/(2*sqrt(Vs)*rm)))], label = "Expected pop size")
		#vline!([Dm/2], label = "Starting deme")
		#plot!([margin]*5, color=:yellow, label = "Deterministic range margin")
		plot(ppf1, grid=false, color=:blue, label="Diploid", legend=false)
		plot!(ppf2, grid=false, color=:green, label="Triploid")
		plot!(ppf3, grid=false, color=:red, label="Tetraploid")
	
		#xlabel!("Space")
		#ylabel!("Population size N")
		
	end every 1
	gif(anim_range, "popsize.gif", fps = 4)
end

# ╔═╡ 822b8dc2-8efb-11eb-3f4d-4b86b1cff707
begin
	function f_trait_agents(h::Habitat)
	trait_agents = []
	cordst = []
	for (i, deme) in enumerate(h.demes)
		for agent in deme.agents
			t = trait(agent)
			p = (i,t)
			push!(cordst,i)
			push!(trait_agents,t)
		end
	end
	trait_agents, cordst
	end
	
	anim_range_trait = @animate for i ∈ 1:300
		sim_habA = sim_hab.data[i]
		trait_means = [trait_mean(deme) for deme in sim_habA.demes]
	    trait_means_p = map(mean, trait_means)
		trait_agents, cordst = f_trait_agents(sim_habA)
		
		plot(g_lin, grid=false, color=:blue, label="Z optimum", linestyle=:dash)
		plot!(cordst,trait_agents, label="Z agents")
		plot!(trait_means_p, grid=false, color=:black, label="Z_mean deme")
		xlabel!("Space")
		ylabel!("Trait Z")
	end every 1
	gif(anim_range_trait, "phenotype.gif", fps = 5)
end

# ╔═╡ Cell order:
# ╠═80bd1070-8e89-11eb-2348-954f56cd5b7b
# ╠═ddff89b0-8e58-11eb-2fef-5b94ca17af89
# ╠═381b957c-e4d7-46bd-b164-c26ce9f3bf7c
# ╠═22341570-8e62-11eb-2dc0-919c8bbaeb55
# ╠═5bc3759e-8f3f-11eb-3458-3d39fe093824
# ╠═5c000970-8f3f-11eb-0c5d-ef59668300a9
# ╠═675ae674-5132-480d-87c5-2e0ad0c6605b
# ╠═8e623f10-8e85-11eb-098f-7184885620f4
# ╠═9f274e80-8e85-11eb-1b80-f5511b4a63f0
# ╠═e34a6b20-8efc-11eb-059d-397e876bd4bc
# ╠═82bcd140-8efb-11eb-1d4a-cf573d3bffe0
# ╠═822b8dc2-8efb-11eb-3f4d-4b86b1cff707
