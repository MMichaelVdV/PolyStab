### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 74d8e4b2-4985-4150-a591-1c58a6907682
using Measures

# ╔═╡ 41fdfc2c-116e-4272-a8b8-ecbdb4f31525
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, DataFrames, StatsPlots

# ╔═╡ 48e0c054-1f38-4517-9afc-558155f24d14
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy, trait_mean, randagent, evolving_neutraldeme

# ╔═╡ 9822f850-ea1a-11eb-2866-9f74ea8160ca
md""" ### Scaling of phenotype and additive genetic variance for different ploidy levels"""

# ╔═╡ 93b286e5-6a06-4026-992c-e1dae2f38b96
md""" We can introduce a parameter `d` to scale the phenotypes for different ploidy levels. Then the genotypic value over all loci can be calculated as follows: $k^{-d}\sum{a}$."""

# ╔═╡ 3265e73d-d2e8-45e1-9f1e-f65c342bd43f
begin
	α=0.5
	L=50
	p=0.5
end

# ╔═╡ fa93b7f4-3bbf-4085-b473-a123c7e6c900
function trait_exp(a::Agent)
	if ploidy(a) == 2
		return (((0.5*(ploidy(a)))^(-a.d))*(sum(a)))
	elseif (ploidy(a) == 4 && a.d == 0.)
		return (((0.5*(ploidy(a)))^(-a.d))*(sum(a)))-(2*(α*L*p))
	elseif (ploidy(a) == 4 && a.d == 0.5)
		return (((0.5*(ploidy(a)))^(-a.d))*(sum(a)))-((2*sqrt(2)-2)*(α*L*p))
	elseif (ploidy(a) == 4 && a.d == 1.)
		return (((0.5*(ploidy(a)))^(-a.d))*(sum(a)))
	end
end	

# ╔═╡ 088dee40-d2e3-4f05-bf5b-26506a406692
#+((sum(a))*((((0.5^a.d)-(0.25^a.d)))*(-0.5*(2-ploidy(a)))))

# ╔═╡ 3fd40211-3c7b-46fd-9040-d41525345558
md""" With this implmentation additive genetic variance ($V_{A}$) scales with a factor $(k/2)^{(1-2d)}$, where `k` is the ploidy level. When we compare the genetic variance of tetraploids to diploids at HWLE with this model for the extreme cases `d`=0 and `d`=1, we get for `d`=0 that the genetic variance is double for tetraploids compared to diploids ($V_{A}$ = $k/2)$ (but note that also the mean phenotypic value shifts). For `d`=1 this is the other way around ($V_{A}$ = $2/k$) (mean phenotypic value is the same of diploids and tetraploids). This is also already evident from the following plots limited to one bi-allelic locus, where the darkness of the dot corresponds to the frequency of the phenotype."""

# ╔═╡ 20ccef7d-6a08-4bbe-ae2f-70c23d7ada83
md"""We can simulate single deme populations for different ploidy levels and different values of `d`. Below I restricted myself to `d`=0, `d`=0.5 and `d`=1."""

# ╔═╡ 3c055079-b69f-4de7-85ea-13283c8d37ac
md""" #### Simulations for d=0"""

# ╔═╡ 18b2bd9f-2702-4c80-93ca-0e637324479b
begin
diploid_d0 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 2500, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 25., Vs = 1.)
tetraploid_d0 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 0., 0., 1.], 200, d=0.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 25., Vs = 1.)
end

# ╔═╡ 45aceacc-aa61-4e2d-a7e3-859866c400c2
begin
Sim_diploid_d0 = evolving_selectiondeme(diploid_d0,trait_exp,1000)#evolving_neutraldeme(diploid_d0,trait_exp,1000) 
Sim_tetraploid_d0 = evolving_selectiondeme(tetraploid_d0,trait_exp,1000)#evolving_neutraldeme(tetraploid_d0,trait_exp,1000) 
end

# ╔═╡ 1927b5f2-fedb-45b7-82b4-2164c0c5e65e
begin
dist0=plot(Normal(mean(Sim_diploid_d0.fta[1]), std(Sim_diploid_d0.fta[1])), color=:black, title="Source population",xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10, linewidth=3,label=false)
#plot!(Normal(mean(Sim_tetraploid_d0.fta[1]), std(Sim_tetraploid_d0.fta[1])), color=:red, label="Tetraploid")
xlabel!("Phenotype")
ylabel!("Frequency")
xlims!(5.,20.)
vline!([12.5],colour="black",linestyle=:solid,label="△z=0",linewidth=3)
vline!([13.5],colour="black",linestyle=:dot,label="△z=1",linewidth=3)
vline!([14.5],colour="black",linestyle=:dash,label="△z=2",linewidth=3)
vline!([15.5],colour="black",linestyle=:dashdot,label="△z=3",linewidth=3)
vline!([16.5],colour="black",linestyle=:dashdotdot,label="△z=4",linewidth=3)
end

# ╔═╡ c4d1336c-ba1e-40a1-a1e4-0354ef020f2f
begin
trait_diploid_d0 = map(mean, Sim_diploid_d0.tm)
trait_tetraploid_d0 = map(mean, Sim_tetraploid_d0.tm)
end

# ╔═╡ 4645f0e8-fa66-4472-a234-b67482ab76c3
begin
	p2d0=plot(trait_diploid_d0, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid, d=0", legend=:bottomright,xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	for (i,t) in enumerate(Sim_diploid_d0.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	plot!(trait_diploid_d0, grid=false, color=:red, label=false,linewidth=3, title="Diploid, d=0", legend=:bottomright,xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([diploid_d0.θ],label=false,colour="black",linestyle=:dash)
	ylims!(5.,50.)
end

# ╔═╡ 3bace483-e385-4257-b057-aa518f8a6859
begin
	p4d0 = plot(trait_tetraploid_d0, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Tetraploid, d=0", legend=:bottomright,xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	for (i,t) in enumerate(Sim_tetraploid_d0.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	plot!(trait_tetraploid_d0, grid=false, color=:red, label=false,linewidth=3, title="Tetraploid, d=0", legend=:bottomright,xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([tetraploid_d0.θ],label=false,colour="black",linestyle=:dash)
	ylims!(5.,50.)
end

# ╔═╡ fda5da73-4a68-4f67-a3d3-47b5570c6973
md""" #### Simulations for d=0.5"""

# ╔═╡ 53f14f64-8b2d-470c-826a-468adb874a81
begin
diploid_d05 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=0.5), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 25., Vs = 1.)
tetraploid_d05 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 0., 0., 1.], 200, d=0.5), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 25., Vs = 1.)
end

# ╔═╡ d83e9e3b-8122-4a01-96d7-38bef868b6d7
trait_exp(diploid_d05.agents[1])

# ╔═╡ 263bb1a6-cb06-4863-b6d0-04963b3eb6bb
trait_exp(tetraploid_d05.agents[1])

# ╔═╡ 50439e5c-e61c-49e1-a373-792f4a9423f4
begin
Sim_diploid_d05 = evolving_selectiondeme(diploid_d05,trait_exp,1000)#evolving_neutraldeme(diploid_d05,trait_exp,1000) 
Sim_tetraploid_d05 = evolving_selectiondeme(tetraploid_d05,trait_exp,1000)#evolving_neutraldeme(tetraploid_d05,trait_exp,1000)
end

# ╔═╡ 223c011f-5f52-4f77-8f17-932fdaf96318
begin
dist05=plot(Normal(mean(Sim_diploid_d05.fta[1]), std(Sim_diploid_d05.fta[1])), color=:blue, label="Diploid", title="d=0.5",xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
plot!(Normal(mean(Sim_tetraploid_d05.fta[1]), std(Sim_tetraploid_d05.fta[1])), color=:red, label="Tetraploid")
xlabel!("Phenotype before selection")
ylabel!("Frequency")
xlims!(15.,35.)
end

# ╔═╡ 0dcc986f-4d06-4c40-8f8e-f5e3f173950b
begin
trait_diploid_d05 = map(mean, Sim_diploid_d05.tm)
trait_tetraploid_d05 = map(mean, Sim_tetraploid_d05.tm)
end

# ╔═╡ 2eba2b3a-3677-4850-b6fd-f88b5ca71d6e
begin
	p2d05=plot(trait_diploid_d05, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid, d=0.5", legend=:bottomright,xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	for (i,t) in enumerate(Sim_diploid_d05.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	plot!(trait_diploid_d05, grid=false, color=:red, label=false,linewidth=3, title="Diploid, d=0.5", legend=:bottomright,xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([diploid_d05.θ],label=false,colour="black",linestyle=:dash)
	ylims!(5.,50.)
end

# ╔═╡ d0c24e83-e1e8-4ec8-ae30-01495bd67440
begin
	p4d05=plot(trait_tetraploid_d05, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Tetraploid, d=0.5", legend=:bottomright,xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	for (i,t) in enumerate(Sim_tetraploid_d05.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	plot!(trait_tetraploid_d05, grid=false, color=:red, label=false,linewidth=3, title="Tetraploid, d=0.5", legend=:bottomright,xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([tetraploid_d05.θ],label=false,colour="black",linestyle=:dash)
	ylims!(5.,50.)
end

# ╔═╡ 767d1302-61ca-44b6-a82a-8ecfdcc0fd62
md""" #### Simulations for d=1"""

# ╔═╡ 9d55b2e6-9335-4894-99b8-6c4fd19e73c6
begin
diploid_d1 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 25., Vs = 1.)
tetraploid_d1 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 0., 0., 1.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 25., Vs = 1.)
end

# ╔═╡ 5919b4da-04ae-439e-a899-65b0f71a440e
begin
Sim_diploid_d1 = evolving_selectiondeme(diploid_d1,trait_exp,1000)#evolving_neutraldeme(diploid_d1,trai
Sim_tetraploid_d1 = evolving_selectiondeme(tetraploid_d1,trait_exp,1000)#evolving_neutraldeme(tetraploid_d1,trait_exp,1000)
end

# ╔═╡ 30489570-43c4-4f8a-bd9d-b7718939d57c
begin
dist1=plot(Normal(mean(Sim_diploid_d1.fta[1]), std(Sim_diploid_d1.fta[1])), color=:blue, label="Diploid", title="d=1",xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
plot!(Normal(mean(Sim_tetraploid_d1.fta[1]), std(Sim_tetraploid_d1.fta[1])), color=:red, label="Tetraploid")
xlabel!("Phenotype before selection")
ylabel!("Frequency")
xlims!(15.,35.)
end

# ╔═╡ 32c2aa17-7ba2-45a8-a69f-5cdde2b1ef0b
begin
trait_diploid_d1 = map(mean, Sim_diploid_d1.tm)
trait_tetraploid_d1 = map(mean, Sim_tetraploid_d1.tm)
end

# ╔═╡ c6f096bc-4ee4-4e26-bc64-660ca910fe42
begin
	p2d1=plot(trait_diploid_d1, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid, d=1", legend=:bottomright,txtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	for (i,t) in enumerate(Sim_diploid_d1.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	plot!(trait_diploid_d1, grid=false, color=:red, label=false,linewidth=3, title="Diploid, d=1", legend=:bottomright,txtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([diploid_d1.θ],label=false,colour="black",linestyle=:dash)
	ylims!(5.,50.)
end

# ╔═╡ 7f2ad115-a7e5-47f7-8c61-afeae195ef32
begin
	p4d1=plot(trait_tetraploid_d1, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Tetraploid, d=1", legend=:bottomright,xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	for (i,t) in enumerate(Sim_tetraploid_d1.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	plot!(trait_tetraploid_d1, grid=false, color=:red, label=false,linewidth=3, title="Tetraploid, d=1", legend=:bottomright,xtickfontsize=10, ytickfontsize=10,xguidefontsize=14, yguidefontsize=14, legendfontsize=10)
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([tetraploid_d1.θ],label=false,colour="black",linestyle=:dash)
	ylims!(5.,50.)
end

# ╔═╡ a3a86212-0766-407e-8be8-b4382b6e1cce
begin
plot(dist0,dist05,dist1,p2d0,p2d05,p2d1,p4d0,p4d05,p4d1,layout=(3,3), size=(1200,800), margins=5mm)
savefig("driftd")
end

# ╔═╡ 87e18c91-b212-4655-9946-b6f62acb361e
begin
	popsize_p2r = map(mean, Sim_diploid_d0.pop)
	popsize_p4r = map(mean, Sim_tetraploid_d0.pop)
	pop0 = plot(popsize_p2r, grid=false, color=:blue, label="Diploid",legend=:bottomright, title="Population growth, d=0",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=16, legendfontsize=10)
	plot!(popsize_p4r, grid=false, color=:red, label="Tetraploid")
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ bec2c142-6531-4159-8d48-4e44a3716cf0
begin
popsize_p2d = map(mean, Sim_diploid_d05.pop)
popsize_p4d = map(mean, Sim_tetraploid_d05.pop)
pop05 = plot(popsize_p2d, grid=false, color=:blue, label="Diploid", title="Population growth, d=0.5",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=16, legendfontsize=10, legend=:bottomright)
plot!(popsize_p4d, grid=false, color=:red, label="Tetraploid")
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ 464e9324-0d67-4a1d-9401-6e4108d1ad9c
begin
popsize_p2a = map(mean, Sim_diploid_d1.pop)
popsize_p4a = map(mean, Sim_tetraploid_d1.pop)
pop1 = plot(popsize_p2a, grid=false, color=:blue, label="Diploid", title="Population growth, d=1",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=16, legendfontsize=10, legend=:bottomright)
plot!(popsize_p4a, grid=false, color=:red, label="Tetraploid")
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ 68c7fcfc-988a-4233-b3f5-89fd71068220
begin
plot(pop0,pop05,pop1, layout=(1,3), size=(1000,300), margins=5mm)
savefig("popd")
end

# ╔═╡ bb053860-bd2e-44d1-842e-238201697628
md""" ##### Cytotype load"""

# ╔═╡ 8764ad8a-ff40-4987-b3ca-5421df3da0b7
begin
	diploid_u0 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 25., Vs = 1.)
	diploid_u05 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 25., Vs = 1.)
	diploid_u10 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.90 0.10 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 25., Vs = 1.)
	diploid_u20 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.80 0.20 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 25., Vs = 1.)
end

# ╔═╡ 37d67078-0f39-4210-9598-ed6a9d6199b2
begin
Sim_diploid_u0 = evolving_selectiondeme(diploid_u0,trait_exp,1000)
Sim_diploid_u05 = evolving_selectiondeme(diploid_u05,trait_exp,1000)
Sim_diploid_u10 = evolving_selectiondeme(diploid_u10,trait_exp,1000)
Sim_diploid_u20 = evolving_selectiondeme(diploid_u10,trait_exp,1000)
end

# ╔═╡ 7414b7d8-950f-48d1-b596-196bf3997f21
begin
popsize_u0 = map(mean, Sim_diploid_u0.pop)
popsize_u05 = map(mean, Sim_diploid_u05.pop)
popsize_010 = map(mean, Sim_diploid_u10.pop)
pop1ug = plot(popsize_u0, grid=false, color=:black, label="u=0", title="Population growth, v=0",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=16, legendfontsize=10, legend=:bottomright)
plot!(popsize_u05, grid=false, color=:red, label="u=0.05")
plot!(popsize_010, grid=false, color=:blue, label="u=0.10")
hline!([diploid_d1.K],label=false,colour="black",linestyle=:dash, linewidth=2)
xlabel!("\$t\$")
ylabel!("Population size")
end

# ╔═╡ cb27b7ec-f9ec-4306-af4a-f9ad2424a82b
begin
	diploid_2u0 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0.90 0. 0.10], θ = 25., Vs = 1.)
	diploid_2u05 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ; 0. 0.95 0. 0.05], θ = 25., Vs = 1.)
	diploid_2u10 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.15 0. 0. ; 0. 0. 0. 0. ; 0. 0.90 0. 0.10], θ = 25., Vs = 1.)
	diploid_2u20 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.70 0.30 0. 0. ; 0. 0. 0. 0. ; 0. 0.70 0. 0.30], θ = 25., Vs = 1.)
end

# ╔═╡ f6b4cbcf-0c46-47b6-859f-4a4012ffa64a
begin
Sim_diploid_2u0 = evolving_selectiondeme(diploid_2u0,trait_exp,250)
#Sim_diploid_2u05 = evolving_selectiondeme(diploid_2u05,trait_exp,1000)
Sim_diploid_2u10 = evolving_selectiondeme(diploid_2u10,trait_exp,250)
#Sim_diploid_2u20 = evolving_selectiondeme(diploid_2u10,trait_exp,1000)
end

# ╔═╡ 7114a8a4-a940-419d-90b3-fd8872fbc00d
begin
popsize_2u0 = map(mean, Sim_diploid_2u0.pop)
#popsize_2u05 = map(mean, Sim_diploid_2u05.pop)
popsize_2010 = map(mean, Sim_diploid_2u10.pop)
#popsize_2020 = map(mean, Sim_diploid_2u20.pop)
pop12ug = plot(popsize_2u0, grid=false, color=:black, label="u=v=0", title="Population growth, u=v",xtickfontsize=10, ytickfontsize=10,xguidefontsize=16,yguidefontsize=16, legendfontsize=10, legend=:bottomright)
#plot!(popsize_2u05, grid=false, color=:red, label="u=v=0.05")
plot!(popsize_2010, grid=false, color=:red, label="u=v=0.15")
#plot!(popsize_2020, grid=false, color=:grey, label="u=v=0.20")
hline!([diploid_d1.K],label=false,colour="black", linestyle=:dash, linewidth=2)
xlabel!("\$t\$")
ylabel!("Population size")
savefig("popcrash.png")
end

# ╔═╡ 926f6e24-2aea-411d-bf70-f1ddc81f04fc
begin
plot(pop1ug,pop12ug, layout=(1,2), size=(1000,300), margins=5mm)
savefig("cytoloadpopsize")
end

# ╔═╡ 1ac4f66d-ef4f-4053-bd46-2e385c14c50a
md""" ##### Functions """

# ╔═╡ 10f262cb-f717-46a3-9c1b-592661b46781
begin
HWE_1(p) = [p,(1-p)]	
HWE_2(p) = [p*p, 2*p*(1-p), (1-p)*(1-p)]
HWE_3(p) = [p*p*p, 3*p*p*(1-p), 3*p*(1-p)*(1-p), (1-p)*(1-p)*(1-p)]
HWE_4(p) = [p*p*p*p, 4*p*p*p*(1-p), 6*p*p*(1-p)*(1-p), 4*p*(1-p)*(1-p)*(1-p), (1-p)*(1-p)*(1-p)*(1-p)]
HWE(p) = [HWE_1(p), HWE_2(p), HWE_3(p), HWE_4(p)]
end

# ╔═╡ 54a62bb2-99d5-4a39-a589-d4306b0d8199
HWEv = HWE(0.5)

# ╔═╡ 729785c1-ba1d-4a21-b2b6-91a0382ba506
genos = [[1.,0.],[2.,1.,0.],[3.,2.,1.,0.],[4.,3.,2.,1.,0.]]

# ╔═╡ 345bdcf6-f9a6-47b9-bb15-2c6d7b99e7ee
begin
	linear_map(v,α,d)=d*α*v.*(2/v[1])
	exponential_map(v,α,d)=α*v.*(2/v[1]).^d
end

# ╔═╡ d3e22239-6f63-4ace-80a4-d6727cab9671
q = reduce(vcat,HWEv)

# ╔═╡ 3f606c7b-c801-46d1-ac3b-64a3cb317c4d
begin
	pHWE0e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="d=0")
	c5 = 1
	for (i,t) in (enumerate(exponential_map.(genos,1.,0.)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c5], ms=7.5)
		c5+= 1
		end
	end
	xlabel!("\$k\$")
	ylabel!("Phenotype")
end

# ╔═╡ d04e8bff-4d74-4dd1-8820-1999f4aeca37
begin
	pHWE25e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="d=0.25")
	c6 = 1
	for (i,t) in (enumerate(exponential_map.(genos,1.,0.25)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c6], ms=7.5)
		c6+= 1
		end
	end
	xlabel!("\$k\$")
	ylabel!("Phenotype")
end

# ╔═╡ f871be81-364e-42e9-9819-f99816e99230
begin
	pHWE50e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="d=0.50")
	c7 = 1
	for (i,t) in (enumerate(exponential_map.(genos,1.,0.50)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c7], ms=7.5)
		c7+= 1
		end
	end
	xlabel!("\$k\$")
	ylabel!("Phenotype")
end

# ╔═╡ 0377f670-822b-44d6-b4c4-6fc804b8fd31
begin
	pHWE75e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="d=0.75")
	c8 = 1
	for (i,t) in (enumerate(exponential_map.(genos,1.,0.75)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c8], ms=7.5)
		c8+= 1
		end
	end
	xlabel!("\$k\$")
	ylabel!("Phenotype")
end

# ╔═╡ 9f39f4db-1629-4727-83e9-08ff0f9f4f4d
begin
	pHWE100e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="d=1")
	c10 = 1
	for (i,t) in (enumerate(exponential_map.(genos,1.,1.)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c10], ms=7.5)
		c10+= 1
		end
	end
	xlabel!("\$k\$")
	ylabel!("Phenotype")
end

# ╔═╡ 9bfb5a50-cea5-4589-b2d2-01d333caed6e
plot(pHWE0e,pHWE25e,pHWE50e,pHWE75e,pHWE100e)

# ╔═╡ Cell order:
# ╟─9822f850-ea1a-11eb-2866-9f74ea8160ca
# ╟─93b286e5-6a06-4026-992c-e1dae2f38b96
# ╠═3265e73d-d2e8-45e1-9f1e-f65c342bd43f
# ╟─fa93b7f4-3bbf-4085-b473-a123c7e6c900
# ╠═088dee40-d2e3-4f05-bf5b-26506a406692
# ╟─3fd40211-3c7b-46fd-9040-d41525345558
# ╠═9bfb5a50-cea5-4589-b2d2-01d333caed6e
# ╟─20ccef7d-6a08-4bbe-ae2f-70c23d7ada83
# ╟─3c055079-b69f-4de7-85ea-13283c8d37ac
# ╠═18b2bd9f-2702-4c80-93ca-0e637324479b
# ╠═d83e9e3b-8122-4a01-96d7-38bef868b6d7
# ╠═263bb1a6-cb06-4863-b6d0-04963b3eb6bb
# ╠═45aceacc-aa61-4e2d-a7e3-859866c400c2
# ╠═1927b5f2-fedb-45b7-82b4-2164c0c5e65e
# ╠═c4d1336c-ba1e-40a1-a1e4-0354ef020f2f
# ╠═4645f0e8-fa66-4472-a234-b67482ab76c3
# ╠═3bace483-e385-4257-b057-aa518f8a6859
# ╟─fda5da73-4a68-4f67-a3d3-47b5570c6973
# ╠═53f14f64-8b2d-470c-826a-468adb874a81
# ╠═50439e5c-e61c-49e1-a373-792f4a9423f4
# ╠═223c011f-5f52-4f77-8f17-932fdaf96318
# ╠═0dcc986f-4d06-4c40-8f8e-f5e3f173950b
# ╠═2eba2b3a-3677-4850-b6fd-f88b5ca71d6e
# ╠═d0c24e83-e1e8-4ec8-ae30-01495bd67440
# ╟─767d1302-61ca-44b6-a82a-8ecfdcc0fd62
# ╠═9d55b2e6-9335-4894-99b8-6c4fd19e73c6
# ╠═5919b4da-04ae-439e-a899-65b0f71a440e
# ╠═30489570-43c4-4f8a-bd9d-b7718939d57c
# ╠═32c2aa17-7ba2-45a8-a69f-5cdde2b1ef0b
# ╠═c6f096bc-4ee4-4e26-bc64-660ca910fe42
# ╠═7f2ad115-a7e5-47f7-8c61-afeae195ef32
# ╠═74d8e4b2-4985-4150-a591-1c58a6907682
# ╠═a3a86212-0766-407e-8be8-b4382b6e1cce
# ╠═87e18c91-b212-4655-9946-b6f62acb361e
# ╠═bec2c142-6531-4159-8d48-4e44a3716cf0
# ╠═464e9324-0d67-4a1d-9401-6e4108d1ad9c
# ╠═68c7fcfc-988a-4233-b3f5-89fd71068220
# ╟─bb053860-bd2e-44d1-842e-238201697628
# ╠═8764ad8a-ff40-4987-b3ca-5421df3da0b7
# ╠═37d67078-0f39-4210-9598-ed6a9d6199b2
# ╠═7414b7d8-950f-48d1-b596-196bf3997f21
# ╠═cb27b7ec-f9ec-4306-af4a-f9ad2424a82b
# ╠═f6b4cbcf-0c46-47b6-859f-4a4012ffa64a
# ╠═7114a8a4-a940-419d-90b3-fd8872fbc00d
# ╠═926f6e24-2aea-411d-bf70-f1ddc81f04fc
# ╟─1ac4f66d-ef4f-4053-bd46-2e385c14c50a
# ╠═41fdfc2c-116e-4272-a8b8-ecbdb4f31525
# ╠═48e0c054-1f38-4517-9afc-558155f24d14
# ╟─10f262cb-f717-46a3-9c1b-592661b46781
# ╟─54a62bb2-99d5-4a39-a589-d4306b0d8199
# ╟─729785c1-ba1d-4a21-b2b6-91a0382ba506
# ╠═345bdcf6-f9a6-47b9-bb15-2c6d7b99e7ee
# ╟─d3e22239-6f63-4ace-80a4-d6727cab9671
# ╟─3f606c7b-c801-46d1-ac3b-64a3cb317c4d
# ╟─d04e8bff-4d74-4dd1-8820-1999f4aeca37
# ╟─f871be81-364e-42e9-9819-f99816e99230
# ╟─0377f670-822b-44d6-b4c4-6fc804b8fd31
# ╟─9f39f4db-1629-4727-83e9-08ff0f9f4f4d
