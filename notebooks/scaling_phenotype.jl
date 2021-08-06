### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 41fdfc2c-116e-4272-a8b8-ecbdb4f31525
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, DataFrames, StatsPlots

# ╔═╡ 48e0c054-1f38-4517-9afc-558155f24d14
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy, trait_mean, randagent

# ╔═╡ 9822f850-ea1a-11eb-2866-9f74ea8160ca
md""" ### Scaling of phenotype and additive genetic variance for different ploidy levels"""

# ╔═╡ 93b286e5-6a06-4026-992c-e1dae2f38b96
md""" We can introduce a parameter `d` to scale the phenotypes for different ploidy levels. Then the genotypic value over all loci can be calculated as follows: $k^{-d}\sum{a}$."""

# ╔═╡ 3265e73d-d2e8-45e1-9f1e-f65c342bd43f
begin
	α=0.1
	L=100
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
diploid_d0 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=0.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.5)
tetraploid_d0 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 0., 0., 1.], 200, d=0.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.5)
end

# ╔═╡ 45aceacc-aa61-4e2d-a7e3-859866c400c2
begin
Sim_diploid_d0 = evolving_selectiondeme(diploid_d0,trait_exp,50)
Sim_tetraploid_d0 = evolving_selectiondeme(tetraploid_d0,trait_exp,50)
end

# ╔═╡ 1927b5f2-fedb-45b7-82b4-2164c0c5e65e
begin
plot(Normal(mean(Sim_diploid_d0.fta[1]), std(Sim_diploid_d0.fta[1])), color=:blue, label="Diploid", title="d=0")
plot!(Normal(mean(Sim_tetraploid_d0.fta[1]), std(Sim_tetraploid_d0.fta[1])), color=:red, label="Tetraploid")
xlabel!("Phenotype before selection")
ylabel!("Frequency")
end

# ╔═╡ c4d1336c-ba1e-40a1-a1e4-0354ef020f2f
begin
trait_diploid_d0 = map(mean, Sim_diploid_d0.tm)
trait_tetraploid_d0 = map(mean, Sim_tetraploid_d0.tm)
end

# ╔═╡ 4645f0e8-fa66-4472-a234-b67482ab76c3
begin
	plot(trait_diploid_d0, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid, d=0", legend=:bottomright)
	for (i,t) in enumerate(Sim_diploid_d0.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([diploid_d0.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ 3bace483-e385-4257-b057-aa518f8a6859
begin
	p2n = plot(trait_tetraploid_d0, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Tetraploid, d=0", legend=:bottomright)
	for (i,t) in enumerate(Sim_tetraploid_d0.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([tetraploid_d0.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ fda5da73-4a68-4f67-a3d3-47b5570c6973
md""" #### Simulations for d=0.5"""

# ╔═╡ 53f14f64-8b2d-470c-826a-468adb874a81
begin
diploid_d05 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=0.5), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = .2)
tetraploid_d05 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 0., 0., 1.], 200, d=0.5), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = .2)
end

# ╔═╡ d83e9e3b-8122-4a01-96d7-38bef868b6d7
trait_exp(diploid_d05.agents[1])

# ╔═╡ 263bb1a6-cb06-4863-b6d0-04963b3eb6bb
trait_exp(tetraploid_d05.agents[1])

# ╔═╡ 50439e5c-e61c-49e1-a373-792f4a9423f4
begin
Sim_diploid_d05 = evolving_selectiondeme(diploid_d05,trait_exp,50)
Sim_tetraploid_d05 = evolving_selectiondeme(tetraploid_d05,trait_exp,50)
end

# ╔═╡ 223c011f-5f52-4f77-8f17-932fdaf96318
begin
plot(Normal(mean(Sim_diploid_d05.fta[1]), std(Sim_diploid_d05.fta[1])), color=:blue, label="Diploid", title="d=0.5")
plot!(Normal(mean(Sim_tetraploid_d05.fta[1]), std(Sim_tetraploid_d05.fta[1])), color=:red, label="Tetraploid")
xlabel!("Phenotype before selection")
ylabel!("Frequency")
end

# ╔═╡ 0dcc986f-4d06-4c40-8f8e-f5e3f173950b
begin
trait_diploid_d05 = map(mean, Sim_diploid_d05.tm)
trait_tetraploid_d05 = map(mean, Sim_tetraploid_d05.tm)
end

# ╔═╡ 2eba2b3a-3677-4850-b6fd-f88b5ca71d6e
begin
	plot(trait_diploid_d05, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid, d=0.5", legend=:bottomright)
	for (i,t) in enumerate(Sim_diploid_d05.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([diploid_d05.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ d0c24e83-e1e8-4ec8-ae30-01495bd67440
begin
	plot(trait_tetraploid_d05, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Tetraploid, d=0.5", legend=:bottomright)
	for (i,t) in enumerate(Sim_tetraploid_d05.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([tetraploid_d05.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ 767d1302-61ca-44b6-a82a-8ecfdcc0fd62
md""" #### Simulations for d=1"""

# ╔═╡ 9d55b2e6-9335-4894-99b8-6c4fd19e73c6
begin
diploid_d1 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.5)
tetraploid_d1 = MixedPloidyDeme(agents = randagent_p(p, α, L, [0., 0., 0., 1.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.5)
end

# ╔═╡ 5919b4da-04ae-439e-a899-65b0f71a440e
begin
Sim_diploid_d1 = evolving_selectiondeme(diploid_d1,trait_exp,50)
Sim_tetraploid_d1 = evolving_selectiondeme(tetraploid_d1,trait_exp,50)
end

# ╔═╡ 30489570-43c4-4f8a-bd9d-b7718939d57c
begin
plot(Normal(mean(Sim_diploid_d1.fta[1]), std(Sim_diploid_d1.fta[1])), color=:blue, label="Diploid", title="d=1")
plot!(Normal(mean(Sim_tetraploid_d1.fta[1]), std(Sim_tetraploid_d1.fta[1])), color=:red, label="Tetraploid")
xlabel!("Phenotype before selection")
ylabel!("Frequency")
end

# ╔═╡ 32c2aa17-7ba2-45a8-a69f-5cdde2b1ef0b
begin
trait_diploid_d1 = map(mean, Sim_diploid_d1.tm)
trait_tetraploid_d1 = map(mean, Sim_tetraploid_d1.tm)
end

# ╔═╡ c6f096bc-4ee4-4e26-bc64-660ca910fe42
begin
	plot(trait_diploid_d1, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid, d=1", legend=:bottomright)
	for (i,t) in enumerate(Sim_diploid_d1.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([diploid_d1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
end

# ╔═╡ 7f2ad115-a7e5-47f7-8c61-afeae195ef32
begin
	plot(trait_tetraploid_d1, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Tetraploid, d=1", legend=:bottomright)
	for (i,t) in enumerate(Sim_tetraploid_d1.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([tetraploid_d1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
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
# ╟─1927b5f2-fedb-45b7-82b4-2164c0c5e65e
# ╠═c4d1336c-ba1e-40a1-a1e4-0354ef020f2f
# ╠═4645f0e8-fa66-4472-a234-b67482ab76c3
# ╠═3bace483-e385-4257-b057-aa518f8a6859
# ╟─fda5da73-4a68-4f67-a3d3-47b5570c6973
# ╠═53f14f64-8b2d-470c-826a-468adb874a81
# ╠═50439e5c-e61c-49e1-a373-792f4a9423f4
# ╟─223c011f-5f52-4f77-8f17-932fdaf96318
# ╠═0dcc986f-4d06-4c40-8f8e-f5e3f173950b
# ╠═2eba2b3a-3677-4850-b6fd-f88b5ca71d6e
# ╠═d0c24e83-e1e8-4ec8-ae30-01495bd67440
# ╟─767d1302-61ca-44b6-a82a-8ecfdcc0fd62
# ╠═9d55b2e6-9335-4894-99b8-6c4fd19e73c6
# ╠═5919b4da-04ae-439e-a899-65b0f71a440e
# ╟─30489570-43c4-4f8a-bd9d-b7718939d57c
# ╠═32c2aa17-7ba2-45a8-a69f-5cdde2b1ef0b
# ╠═c6f096bc-4ee4-4e26-bc64-660ca910fe42
# ╠═7f2ad115-a7e5-47f7-8c61-afeae195ef32
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
