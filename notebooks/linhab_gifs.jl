### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 80bd1070-8e89-11eb-2348-954f56cd5b7b
using Statistics, Plots 

# ╔═╡ ddff89b0-8e58-11eb-2fef-5b94ca17af89
using PolyStab: AbstractDeme, randagent_p, MixedPloidyDeme, Habitat, linear_gradient, initiate_habitat, evolving_habitat, ploidy_freq, heterozygosities_p, trait_mean, trait

# ╔═╡ 3a757fed-7240-41f9-8533-b62490d54d4c
using StatsBase, StatsPlots

# ╔═╡ a14f04bc-a720-49ea-b2b5-421da6412a52
using DataFrames, GLM

# ╔═╡ 3ebcf7ed-777a-4631-944d-e5a839302f23
using Measures

# ╔═╡ 2f0e9dcb-c6b1-43c3-8512-736e92f45ece
md""" ##### Simulations with mixed ploidy populations along a linear gradient."""

# ╔═╡ 381b957c-e4d7-46bd-b164-c26ce9f3bf7c
d = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.], 0), K=25,  OV = [1.0 0. 0.0 0.0; 0. 1.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], UG = [0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0; 0. 0. 0. 0.; 0. 1. 0.0 0.], Vs=1.)

# ╔═╡ 22341570-8e62-11eb-2dc0-919c8bbaeb55
habi = Habitat(demes=[d], Dm=250., b=.1, σ=0.5)

# ╔═╡ 5bc3759e-8f3f-11eb-3458-3d39fe093824
g_lin = linear_gradient(.1,12.5,250)

# ╔═╡ 5c000970-8f3f-11eb-0c5d-ef59668300a9
hab = initiate_habitat(d, g_lin, 0.5, 0.5, 50, 25, [0., 1., 0., 0.])

# ╔═╡ 675ae674-5132-480d-87c5-2e0ad0c6605b
habs = Habitat(demes=hab.demes, σ=0.5, b=.1, θ=12.5, Dm=250.)

# ╔═╡ e72ca08d-7c01-43a4-8b23-8adfe822f671
begin
d1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.], 0), K=25,  OV = [1.0 0.1 0.0 0.0; 0.1 1.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], UG = [0.0 0.0 0.0 0.0; 0.95 0.05 0.0 0.0; 0. 0. 0. 0.; 0. 1. 0.0 0.], Vs=1.2, θ=10.)
d2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.], 25), K=25,  OV = [1.0 0.1 0.0 0.0; 0.1 1.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], UG = [0.0 0.0 0.0 0.0; 0.95 0.05 0.0 0.0; 0. 0. 0. 0.; 0. 1. 0.0 0.], Vs=1.2, θ=12.5)
d3 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.], 0), K=25,  OV = [1.0 0.1 0.0 0.0; 0.1 1.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], UG = [0.0 0.0 0.0 0.0; 0.95 0.05 0.0 0.0; 0. 0. 0. 0.; 0. 1. 0.0 0.], Vs=1.2, θ=15.)
end

# ╔═╡ afa9029a-53cf-46b5-9ab5-b9bffe99ef47
#habs = Habitat(demes=[d1,d2,d3])

# ╔═╡ 8e623f10-8e85-11eb-098f-7184885620f4
sim_hab = evolving_habitat(habs, 500)

# ╔═╡ 10f9f90b-df70-4a6d-afcc-6b2ee4fdf68d
sim_start = evolving_habitat(habs, 0)

# ╔═╡ 9f274e80-8e85-11eb-1b80-f5511b4a63f0
begin
	K = d.K
	σ = hab.σ
	b = hab.b
	Vs = d.Vs
	rm = d.rm
	Dm = hab.Dm
	s = 1
	pop_sizes_start = [length(deme) for deme  in sim_start[1].demes]
	pop_sizes = [length(deme) for deme  in sim_hab[1].demes]
	margin = (sqrt(2)*b*σ)/((2*rm*sqrt(Vs))-b*σ) .>= 0.15.*pop_sizes.*σ*sqrt(s)
	ppf1 = [ploidy_freq(deme)[2] for deme  in sim_hab[1].demes]
	ppf2 = [ploidy_freq(deme)[3] for deme  in sim_hab[1].demes]
	ppf3 = [ploidy_freq(deme)[4] for deme  in sim_hab[1].demes]
	pop_sizes = map(Statistics.mean, pop_sizes)
	p1 = plot(pop_sizes, grid=false, color=:black, label=false,xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=16, legendfontsize=10, linewidth=2, title="Diploid")
	
	plot!(pop_sizes, grid=false, color=:black, label=false, legendfontsize=5,legend=false)
	plot!(pop_sizes_start, grid=false, color=:black, label=false, legendfontsize=5,linestyle=:dash,legend=false,linewidth=2)
	hline!([K], label = "K", color=:black,linestyle=:dash, linewidth=2)
	#hline!([K*(1-(σ*b)*(1/(2*sqrt(Vs)*rm)))], label = "Expected pop size")
	#vline!([Dm/2], label = "Starting deme", color=:black,linestyle=:dash)
	#plot!([margin]*10, label = "Deterministic range margin")
	#plot!(ppf1, grid=false, color=:blue, label="Diploid",linewidth=2)
	#plot!(ppf2, grid=false, color=:purple, label="Triploid")
	#plot!(ppf3, grid=false, color=:red, label="Tetraploid",linewidth=2)
	xlabel!("Space")
	ylabel!("Population size")
	#savefig(p1, "popsize1D")
end

# ╔═╡ ca5f656c-1c71-4fe8-84e8-085ea64338ed
begin
	popsize(d::AbstractDeme) = length(d.agents)

	get_first(v)=v[1]
	get_second(v)=v[2]
	get_third(v)=v[3]
	get_fourth(v)=v[4]

	pop_sizesb = popsize.(sim_hab[1].demes)
	ppf1b = get_second.(ploidy_freq.(sim_hab[1].demes))
	ppf2b = get_third.(ploidy_freq.(sim_hab[1].demes))
	ppf3b = get_fourth.(ploidy_freq.(sim_hab[1].demes))
	ppf4b = argmax.(ploidy_freq.(sim_hab[1].demes))
	
	ph2d = heatmap(reshape(pop_sizesb, 1, length(pop_sizesb)), c=cgrad([:white,:grey,:black]))
	p1h2d = heatmap(reshape(ppf1b, 1, length(ppf1b)), c=cgrad([:white,:blue]))
	p2h2d = heatmap(reshape(ppf2b, 1, length(ppf2b)), c=cgrad([:white,:green]))
	p3h2d = heatmap(reshape(ppf3b, 1, length(ppf3b)), c=cgrad([:white,:red]))
	p4h2d = heatmap(reshape(ppf4b, 1, length(ppf4b)), c=cgrad([:white,:blue,:green,:red]), clim= (0.,4.))

end

# ╔═╡ 30ebd7fd-2c6b-4e53-a664-be146460bc90
plot(ph2d,p1h2d,p2h2d,p3h2d,p4h2d, layout = (5,1))

# ╔═╡ 24e1603f-8934-4a8f-a5c5-216f5dedef46
[(popsize.(sim_hab[1].demes))]

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

# ╔═╡ 1aae568d-9e9a-48c5-ab00-47f989abc1fa
begin
	het_demes, cordsh = f_het_demes(sim_hab[1])
	het_demes_start, cordsh_start = f_het_demes(sim_start[1])
		
	p1h = plot(cordsh, het_demes, grid=false, color=:black, label="Vg_mean deme",linewidth=2,xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=16, legendfontsize=10,legend=false)
	p1h_start = plot!(cordsh, het_demes_start, grid=false, color=:black, label="Vg_mean deme",linestyle=:dash,linewidth=2)
	#vline!([Dm/2], label = "Starting deme")
		#plot!([margin]*1, label = "Deterministic range margin")
		#plot!(ppf1, grid=false, color=:blue, label="Diploid")
		#plot!(ppf2, grid=false, color=:green, label="Triploid")
		#plot!(ppf3, grid=false, color=:red, label="Tetraploid")
		
		
	xlabel!("Space")
	ylabel!("Genetic variance")
	ylims!(0.,0.7)

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

# ╔═╡ faf81b0e-c7d4-4b54-a74a-143c63443089
begin

		trait_means = [trait_mean(deme) for deme in sim_hab[1].demes]
		trait_means_start = [trait_mean(deme) for deme in sim_start[1].demes]
	
	    trait_means_p = map(mean, trait_means)
		trait_agents, cordst = f_trait_agents(sim_hab[1])
	
		trait_means_p_start = map(mean, trait_means_start)
		trait_agents_start, cordst_start = f_trait_agents(sim_start[1])
	
		
		p1t = plot(g_lin, grid=false, color=:black, label="Z optimum", linestyle=:dash,linewidth=2,xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=16, legendfontsize=10,legend=false)
		plot!(cordst,trait_agents, label="Z agents", color=:black)
		plot!(trait_means_p, grid=false, color=:black, label="Z_mean deme",linewidth=2)
		plot!(cordst_start,trait_agents_start, label="Z agents", color=:blue, linewidth=2)
		xlabel!("Space")
		ylabel!("Phenotype")

end

# ╔═╡ 5df4d750-d4b9-4d05-9a94-0eba9b3e1165
begin
plot(p1,p1h,p1t,layout=(3,1), size=(400,1000), margin=5mm)
savefig("linhabdiploid.png")
end

# ╔═╡ 598b6c2c-9eeb-4907-932e-f359839b06f5
md""" ##### Polyploid estabishment simulations """

# ╔═╡ ef0908aa-a2a5-46ba-af69-d60ed9738ed7
function grid_search(t)
	ploidy = []
	param = []
	pop_size = []

	for bs in range(0.1, stop=1., length=t)
		for rep in 1:10
		d = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.], 0, d=1.), K=25,  OV = [1.0 0. 0.0 0.0; 0. 1.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], UG = [0.0 0.0 0.0 0.0; 0.95 0.05 0.0 0.0; 0. 0. 0. 0.; 0. 1. 0.0 0.], Vs=0.5)
		habi = Habitat(demes=[d], Dm=250., b=bs, σ=0.5)
		g_lin = linear_gradient(bs,12.5,250)
		hab = initiate_habitat(d, g_lin, 0.5, 0.5, 50, 25, [0., 1., 0., 0.])
		habs = Habitat(demes=hab.demes, σ=0.5, b=bs, θ=12.5, Dm=250.)
		sim_hab = evolving_habitat(habs, 250)
			
		#frequency of 4n > 2n demes
		pop_sizesb = popsize.(sim_hab[1].demes) # > 0
		ppf2b = get_second.(ploidy_freq.(sim_hab[1].demes))
		ppf3b = get_third.(ploidy_freq.(sim_hab[1].demes))
		ppf4b = get_fourth.(ploidy_freq.(sim_hab[1].demes))
		ppfdb = argmax.(ploidy_freq.(sim_hab[1].demes)) # 2 or 4
			
		
		pestab = sum(popsize.(sim_hab[1].demes)) .> 0 ? sum((argmax.(ploidy_freq.(sim_hab[1].demes)) .== 4))/sum(popsize.(sim_hab[1].demes) .> 0) : 0.
		
		push!(ploidy,pestab)
		push!(param, bs)
		pop = sum(popsize.(sim_hab[1].demes) .> 0) #number of populated demes
		push!(pop_size, pop)
		end
	end
	ploidy, param, pop_size
end

# ╔═╡ f2f02669-7701-4c0b-8336-2b35016b8d58
g = [i for i in range(0.1, stop=1., length=10)]

# ╔═╡ b78f867a-8aac-47da-aed6-5598f641228c
sum(popsize.(sim_hab[1].demes) .> 0)

# ╔═╡ 34af6157-5873-4eff-9a6b-53b803f08c63
sum((argmax.(ploidy_freq.(sim_hab[1].demes)) .== 4))/sum((argmax.(ploidy_freq.(sim_hab[1].demes)) .== 2))

# ╔═╡ 617c281b-7f5f-4ba1-8420-e99a80ef6be1
sim_b = grid_search(10)

# ╔═╡ e4bedab0-051e-41ae-bd09-291e09ec826a
sim_b[1]

# ╔═╡ 7bf93c00-7f82-4d81-a469-4624dd971583
[(sim_b[1][i:i+9]) for i in 1:10:91]*1.

# ╔═╡ 6a1e921e-bf21-4781-b502-9659748f2b3b
[(sim_b[3][i:i+9]) for i in 1:10:91]*1.

# ╔═╡ dcc3e98d-0d5c-4009-b425-5dd055c3ba0e
begin
binned_bm = [mean(sim_b[1][i:i+9]) for i in 1:10:91]
binned_be = [sem(sim_b[1][i:i+9]) for i in 1:10:91]
end

# ╔═╡ 3c3afc15-33e6-44f8-80d6-a533b03c0b10
binned_bm

# ╔═╡ 44d955e1-e475-4c6e-9e33-9de9fdd26f66


# ╔═╡ d9d3b3fd-00d7-4b37-a9b6-dd9a378c55b8


# ╔═╡ 1bdadef7-19eb-4e9d-8e0c-727f62ff00ca
sim_b[3]

# ╔═╡ 48a184dc-e285-49af-9026-0cd91ee420a4
begin
	binned_dm = [mean(sim_b[3][i:i+9]) for i in 1:10:91]
	binned_de = [sem(sim_b[3][i:i+9]) for i in 1:10:91]
end

# ╔═╡ 54cc6552-910d-4eb5-adf3-e2d3f0a9ce49
begin
SE=g
crit_SE=binned_dm           
df_SE=DataFrame([SE crit_SE])
rename!(df_SE,:x1 => :SE)
rename!(df_SE,:x2 => :crit_SE)
end

# ╔═╡ e6f133d3-8051-4196-8203-bab6d42fffd8
fit_SE = @formula(crit_SE ~ SE)

# ╔═╡ feb13c26-08c4-48ce-b0ca-d3df559fdadd
reg_SE = lm(fit_SE, df_SE)

# ╔═╡ 9873021d-6a2c-48ad-b4a9-35690495eaec
r2(reg_SE)

# ╔═╡ b073cffa-a256-4c0e-aa06-559e7e43b94d
begin
t1 = scatter(g,binned_bm,yerror=binned_be,colour=:black,title=:"u=0.05, σ=0.5",label=:false, xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=14, legendfontsize=12)
xlabel!("b")
ylabel!("P estab")
ylims!((0.0,1.))
end

# ╔═╡ 3f581fc3-41a7-4a7e-945a-048cbc637d55
begin
t2 = scatter(g,binned_dm,yerror=binned_de,colour=:black,title=:"u=0.05, σ=0.5",label=:false, xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=14, legendfontsize=14)
xlabel!("b")
ylabel!("Number of populated demes")
ylims!((0.0,80.))
	
lr2(x) = 61.7667 -54.0303*x
plot!([0.:0.005:1. ...],lr2.([0.:0.005:1. ...]), colour =:black, linewidth=2, label="\$y=91.147-82.849x, R^2=0.969\$")
	
end

# ╔═╡ a7b913f6-9342-4cd6-bbce-1ae9bbfa1922
begin
	plot(t1,t2, size=(1200,500), margin=5mm)
	#savefig("linhabp2")
end

# ╔═╡ Cell order:
# ╟─2f0e9dcb-c6b1-43c3-8512-736e92f45ece
# ╠═80bd1070-8e89-11eb-2348-954f56cd5b7b
# ╠═ddff89b0-8e58-11eb-2fef-5b94ca17af89
# ╠═381b957c-e4d7-46bd-b164-c26ce9f3bf7c
# ╠═22341570-8e62-11eb-2dc0-919c8bbaeb55
# ╠═5bc3759e-8f3f-11eb-3458-3d39fe093824
# ╠═5c000970-8f3f-11eb-0c5d-ef59668300a9
# ╠═675ae674-5132-480d-87c5-2e0ad0c6605b
# ╠═e72ca08d-7c01-43a4-8b23-8adfe822f671
# ╠═afa9029a-53cf-46b5-9ab5-b9bffe99ef47
# ╠═8e623f10-8e85-11eb-098f-7184885620f4
# ╠═10f9f90b-df70-4a6d-afcc-6b2ee4fdf68d
# ╠═9f274e80-8e85-11eb-1b80-f5511b4a63f0
# ╠═1aae568d-9e9a-48c5-ab00-47f989abc1fa
# ╠═faf81b0e-c7d4-4b54-a74a-143c63443089
# ╠═5df4d750-d4b9-4d05-9a94-0eba9b3e1165
# ╠═e34a6b20-8efc-11eb-059d-397e876bd4bc
# ╠═ca5f656c-1c71-4fe8-84e8-085ea64338ed
# ╠═30ebd7fd-2c6b-4e53-a664-be146460bc90
# ╠═24e1603f-8934-4a8f-a5c5-216f5dedef46
# ╠═82bcd140-8efb-11eb-1d4a-cf573d3bffe0
# ╠═822b8dc2-8efb-11eb-3f4d-4b86b1cff707
# ╟─598b6c2c-9eeb-4907-932e-f359839b06f5
# ╠═ef0908aa-a2a5-46ba-af69-d60ed9738ed7
# ╠═f2f02669-7701-4c0b-8336-2b35016b8d58
# ╠═b78f867a-8aac-47da-aed6-5598f641228c
# ╠═34af6157-5873-4eff-9a6b-53b803f08c63
# ╠═617c281b-7f5f-4ba1-8420-e99a80ef6be1
# ╠═e4bedab0-051e-41ae-bd09-291e09ec826a
# ╠═3a757fed-7240-41f9-8533-b62490d54d4c
# ╠═7bf93c00-7f82-4d81-a469-4624dd971583
# ╠═6a1e921e-bf21-4781-b502-9659748f2b3b
# ╠═dcc3e98d-0d5c-4009-b425-5dd055c3ba0e
# ╠═3c3afc15-33e6-44f8-80d6-a533b03c0b10
# ╠═44d955e1-e475-4c6e-9e33-9de9fdd26f66
# ╠═d9d3b3fd-00d7-4b37-a9b6-dd9a378c55b8
# ╠═1bdadef7-19eb-4e9d-8e0c-727f62ff00ca
# ╠═48a184dc-e285-49af-9026-0cd91ee420a4
# ╠═a14f04bc-a720-49ea-b2b5-421da6412a52
# ╠═54cc6552-910d-4eb5-adf3-e2d3f0a9ce49
# ╠═e6f133d3-8051-4196-8203-bab6d42fffd8
# ╠═feb13c26-08c4-48ce-b0ca-d3df559fdadd
# ╠═9873021d-6a2c-48ad-b4a9-35690495eaec
# ╠═b073cffa-a256-4c0e-aa06-559e7e43b94d
# ╠═3f581fc3-41a7-4a7e-945a-048cbc637d55
# ╠═3ebcf7ed-777a-4631-944d-e5a839302f23
# ╠═a7b913f6-9342-4cd6-bbce-1ae9bbfa1922
