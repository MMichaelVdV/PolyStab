### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 8ea08d1e-8bff-11eb-3fbb-f74cf53a2eed
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, DataFrames, StatsPlots

# ╔═╡ ded34c80-8699-11eb-33ce-214cb4f59699
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, evolving_neutraldeme, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy, trait_mean, randagent, evolving_haploiddeme, mate_p, evolving_selectiondemeh

# ╔═╡ 91826552-6b17-491f-ad1b-77e06dfcdeb6
using PolyStab: malthusian_fitness, trait_add, trait_dom, trait_rec

# ╔═╡ f4182f59-e44c-4f1d-838b-f4cd9ae8f22d
using PolyStab: ploidy_freq, f_trait_agents, mating_PnB_x, mutate

# ╔═╡ b79eb4e0-8696-11eb-0949-555c0c3c411f
md"""### Single deme dynamics of a mixed ploidy population"""

# ╔═╡ 09e48960-9875-11eb-35fe-5122ae8dbe5f
md""" ##### Initialization of demes with different ploidy levels (haploid, diploid, tetraploid)"""

# ╔═╡ 217c5170-8e52-11eb-175c-c9b8476ae81b
d_p1 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [1., 0., 0., 0.],150), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 0.1)

# ╔═╡ 3f70dcc0-8bfd-11eb-3df6-6709900e71ce
d_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [0., 1., 0., 0.],150), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 0.1)

# ╔═╡ 39b5e750-8dbd-11eb-3f64-c12fa3fe01d5
d_p4 = MixedPloidyDeme(agents = randagent_p(0.5, 0.2, 200, [0., 0., 0., 1.],150), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 0.1)

# ╔═╡ 3106b1d0-9875-11eb-36ac-0d2b98a73fd4
md""" ##### Dynamics in a neutral deme (random mating, no selection)"""

# ╔═╡ a8d89200-8e52-11eb-1277-b1857f7ea0fe
begin
neutral_p1 = evolving_haploiddeme(d_p1, 1000)
neutral_p2 = evolving_neutraldeme(d_p2, 1000)
neutral_p4 = evolving_neutraldeme(d_p4, 1000)
end

# ╔═╡ 66a25c90-8bff-11eb-2a53-09b98728f17c
begin
	Hₒ_p1 = map(mean, neutral_p1.het)
	Hₒ_p2 = map(mean, neutral_p2.het)
	Hₒ_p4 = map(mean, neutral_p4.het)
	expected_heterozygosity(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	eh_1(H₀, t, N) = ((1.0-1.0/(N))^t)*H₀
	eh_2(H₀, t, N) = ((1.0-1.0/(2*N))^t)*H₀
	eh_4(H₀, t, N) = ((1.0-1.0/(4*N))^t)*H₀
	plot(Hₒ_p1, grid=false, color=:black, label="\$H_op1(t)\$", title="LOH")
	plot!(Hₒ_p2, grid=false, color=:blue, label="\$H_op2(t)\$")
	plot!(Hₒ_p4, grid=false, color=:red, label="\$H_op4(t)\$")
	#plot!(0:sim_p2.ngen+1, 
		#t->expected_heterozygosity(0.5*(1-0.5), t, 50),
		#linestyle=:dash, color=:black, 
		#label = "\$(1-1/N)^t H_o(0)\$")
	plot!(0:neutral_p1.ngen+1, 
		t->eh_1(0.5*(1-0.5), t, 150),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	plot!(0:neutral_p2.ngen+1, 
		t->eh_2(0.5*(1-0.5), t, 150),
		linestyle=:dash, color=:blue, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	plot!(0:neutral_p4.ngen+1, 
		t->eh_4(0.5*(1-0.5), t, 150),
		linestyle=:dash, color=:red, 
		label = "\$(1-1/4N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 621f41a0-9885-11eb-2bf8-a530b5cc7e86
begin
n = 200 #number of loci
k = 2 #ploidy
#α = 1/√(k*n) 
α = .2 
a2 = randagent(0.5, α, n, 2)
b2 = randagent(0.5, α, n, 2)
o2 = map(x->trait(mate_p(a2,b2)), 1:2500) 

distro = map(x->trait(randagent(0.5, α, n, 2)),1:2500)

O_mean2 = mean(o2)
E_mean2 = 0.5*(trait(a2) + trait(b2)) 
O_var2 = var(o2)
E_var2 = 0.5 * n * (α^2 / 4)
	
a4 = randagent(0.5, α, n, 4)
b4 = randagent(0.5, α, n, 4)
o4 = map(x->trait(mate_p(a4,b4)), 1:2500) 

O_mean4 = mean(o4)
E_mean4 = 0.5*(trait(a4) + trait(b4)) 
O_var4 = var(o4)
E_var4 = 
	
plot(Normal(O_mean2, sqrt(O_var2)),color=:blue, label="Observed 2N")
plot!(Normal(mean(distro), std(distro)), label="HWLE 2N")
plot!(Normal(O_mean4, sqrt(O_var4)),color=:red, label="Observed 4N")
plot!(Normal(E_mean2,sqrt(E_var2)),color=:black, label="Expected 2N")
xlabel!("\$phenotype\$")
ylabel!("\$\$")
end

# ╔═╡ a6806a20-9927-11eb-3fa6-3769de70a1f5
begin
	traitneutral_p1 = map(mean, neutral_p1.tm)
	traitneutral_p2 = map(mean, neutral_p2.tm)
	traitneutral_p4 = map(mean, neutral_p4.tm)
	plot(traitneutral_p1, grid=false, color=:black, label=false)
	plot!(traitneutral_p2, grid=false, color=:blue, label=false)
	plot!(traitneutral_p4, grid=false, color=:red, label=false)
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 3ec21690-992c-11eb-3b9f-4143af21edc5
begin
	p1n = plot(traitneutral_p1, grid=false, color=:red,linewidth=3,label="Mean phenotype", title="Haploid", legend=:bottomright)
	for (i,t) in enumerate(neutral_p1.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(17,25)
end

# ╔═╡ 415089a2-992c-11eb-24c9-ff5fa01c0273
begin
	p2n = plot(traitneutral_p2, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Diploid", legend=:bottomright)
	for (i,t) in enumerate(neutral_p2.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(17,25)
end

# ╔═╡ 42f31610-992c-11eb-03e1-5b5f8ca6a2fc
begin
	p4n = plot(traitneutral_p4, grid=false, color=:red, label="Mean phenotype",linewidth=3, title="Tetraploid",legend=:bottomright)
	for (i,t) in enumerate(neutral_p4.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Phenotype")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(17,25)
end

# ╔═╡ 1b1541de-8dbf-11eb-041c-815c116efeeb
begin
	pop_p1 = map(mean, neutral_p1.p1)
	pop_p2 = map(mean, neutral_p2.p2)
	pop_p4 = map(mean, neutral_p4.p4)
	plot(pop_p1, grid=false, color=:black, label="\$pop_p1(t)\$", title="Population size")
	plot!(pop_p2, grid=false, color=:blue, label="\$pop_p2(t)\$")
	plot!(pop_p4, grid=false, color=:red, label="\$pop_p4(t)\$")

	xlabel!("\$t\$")
	ylims!(145,155)
	ylabel!("\$pop(t)\$")
end

# ╔═╡ 3fc09382-986b-11eb-29d6-57eeb0eddb77
md""" ##### Dynamics in a deme with densitiy dependence and stabilizing selection"""

# ╔═╡ 47dc60c0-9876-11eb-3d70-cd03509cbad8
begin
stabsel_p1 = evolving_selectiondemeh(d_p1,500)
stabsel_p2 = evolving_selectiondeme(d_p2,500)
stabsel_p4 = evolving_selectiondeme(d_p4,500)
end

# ╔═╡ d6c03b40-9939-11eb-244b-3d747272f440
begin
plot(Normal(mean(stabsel_p1.fta[1]), std(stabsel_p1.fta[1])), color=:black, label="Haploid")
plot!(Normal(mean(stabsel_p2.fta[1]), std(stabsel_p2.fta[1])), color=:blue, label="Diploid")
plot!(Normal(mean(stabsel_p4.fta[1]), std(stabsel_p4.fta[1])), color=:red, label="Tetraploid")
xlabel!("Phenotype before selection")
end

# ╔═╡ fc96ac10-9884-11eb-1735-d751fd99c805
begin
	Hₒs_p1 = map(mean, stabsel_p1.het)
	Hₒs_p2 = map(mean, stabsel_p2.het)
	Hₒs_p4 = map(mean, stabsel_p4.het)

	plot(Hₒs_p1, grid=false, color=:black, label="\$H_op1(t)\$", title="LOH")
	plot!(Hₒs_p2, grid=false, color=:blue, label="\$H_op2(t)\$")
	plot!(Hₒs_p4, grid=false, color=:red, label="\$H_op4(t)\$")
	#plot!(0:sim_p2.ngen+1, 
		#t->expected_heterozygosity(0.5*(1-0.5), t, 50),
		#linestyle=:dash, color=:black, 
		#label = "\$(1-1/N)^t H_o(0)\$")
	plot!(0:neutral_p1.ngen+1, 
		t->eh_1(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	plot!(0:neutral_p2.ngen+1, 
		t->eh_2(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:blue, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	plot!(0:neutral_p4.ngen+1, 
		t->eh_4(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:red, 
		label = "\$(1-1/4N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 66257e10-9879-11eb-2367-d706689f326a
begin
	traitmean_p1 = map(mean, stabsel_p1.tm)
	traitmean_p2 = map(mean, stabsel_p2.tm)
	traitmean_p4 = map(mean, stabsel_p4.tm)
	plot(traitmean_p1, grid=false, color=:black, label=false)
	plot!(traitmean_p2, grid=false, color=:blue, label=false)
	plot!(traitmean_p4, grid=false, color=:red, label=false)
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 0501879e-9928-11eb-01c1-23fdad736664
begin
	p1 = plot(traitmean_p1, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright, title="Haploid")
	for (i,t) in enumerate(stabsel_p1.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(18,22)
end

# ╔═╡ 70a06c90-992a-11eb-1a39-ff907051a107
begin
	p2 = plot(traitmean_p2, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright, title="Diploid")
	for (i,t) in enumerate(stabsel_p2.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(18,22)
end

# ╔═╡ 7209e7a0-992a-11eb-0690-159af85e9c8d
begin
	p4 = plot(traitmean_p4, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright, title="Tetraploid")
	for (i,t) in enumerate(stabsel_p4.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(18,22)
end

# ╔═╡ 63ad1c50-987f-11eb-1512-df6fc858aaca
begin
traitp1 = plot(traitmean_p1, grid=false, color=:black, label=false)
traitp2 = plot(traitmean_p2, grid=false, color=:blue, label=false)
traitp3 = plot(traitmean_p4, grid=false, color=:red, label=false)
plot(traitp1,traitp2,traitp3, layout = (3, 1))
end

# ╔═╡ 44614280-9920-11eb-302b-6f959c8c2b94
begin
	plot(Normal(mean(stabsel_p1.fta[1]), std(stabsel_p1.fta[1])), label="", title="Haploid")
	plot!(Normal(mean(stabsel_p1.fta[end]), std(stabsel_p1.fta[end])), label="")
	xlabel!("Phenotype")
end

# ╔═╡ 809cf8b2-9926-11eb-2b5e-d97d23d44531
begin
	histogram(stabsel_p1.fta[1], bins = 25, fillalpha = 0.4, title="Haploid start")
	xlabel!("Phenotype")
end

# ╔═╡ 651f5aa2-9927-11eb-392f-c7944c0ad0aa
begin
	histogram(stabsel_p1.fta[end], bins = 25, fillalpha = 0.4, title="Haploid after")
	xlabel!("Phenotype")
end

# ╔═╡ 134a9ebe-9921-11eb-3d7e-cdd5d3ff9ee1
begin
	plot(Normal(mean(stabsel_p2.fta[1]), std(stabsel_p2.fta[1])), label="", title="Diploid")
	plot!(Normal(mean(stabsel_p2.fta[end]), std(stabsel_p2.fta[end])), label="")
	xlabel!("Phenotype")
end

# ╔═╡ bae32f30-9926-11eb-3f23-cf639e204a0e
begin
	histogram(stabsel_p2.fta[1], bins = 25, fillalpha = 0.4, title="Diploid start")
	xlabel!("Phenotype")
end

# ╔═╡ c1e57720-9926-11eb-1d2d-058c7c88e98d
begin
	histogram(stabsel_p2.fta[end], bins = 25, fillalpha = 0.4, title="Diploid after")
	xlabel!("Phenotype")
end

# ╔═╡ 158facc0-9921-11eb-25ef-a759b775ecdd
begin
	plot(Normal(mean(stabsel_p4.fta[1]), std(stabsel_p4.fta[1])), label="",title="Tetraploid")
	plot!(Normal(mean(stabsel_p4.fta[end]), std(stabsel_p4.fta[end])), label="")
	xlabel!("Phenotype")
end

# ╔═╡ cb25c830-9926-11eb-01f3-b3ddd767318b
begin
	histogram(stabsel_p4.fta[1], bins = 25, fillalpha = 0.4, title="Tetraploid start")
	xlabel!("Phenotype")
end

# ╔═╡ 32ea0626-9019-4f27-8844-70c4f696549e
begin
	histogram(stabsel_p4.fta[end], bins = 25, fillalpha = 0.4, title="Tetraploid after")
	xlabel!("Phenotype")
end

# ╔═╡ db982b20-9923-11eb-285b-b9d57489b052
md""" Phenotypic variance seems to persist longer in tetraploids under stabilizing selection, this coulbd be mainly due to the effect of less genetic drift or something else. Tetraploid population takes longer to reach equilibrium."""

# ╔═╡ c9b2b8d0-9879-11eb-294a-952184f8bf9c
begin
	popsize_p1 = map(mean, stabsel_p1.pop)
	popsize_p2 = map(mean, stabsel_p2.pop)
	popsize_p4 = map(mean, stabsel_p4.pop)
	plot(popsize_p1, grid=false, color=:black, label=false)
	plot!(popsize_p2, grid=false, color=:blue, label=false)
	plot!(popsize_p4, grid=false, color=:red, label=false)
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ 9cdac49e-9884-11eb-2be7-832739ca5d66
begin
	popp1 = plot(popsize_p1, grid=false, color=:black, label=false)
	popp2 = plot(popsize_p2, grid=false, color=:blue, label=false)
	popp4 = plot(popsize_p4, grid=false, color=:red, label=false)
	plot(popp1,popp2,popp4, layout = (3, 1))
end

# ╔═╡ fe5aa0f0-9939-11eb-1416-91041276c925
md""" ###### Genetic variance"""

# ╔═╡ 5ca18050-9955-11eb-2703-47174bc221a8
begin
HWE_1(p) = [p,(1-p)]	
HWE_2(p) = [p*p, 2*p*(1-p), (1-p)*(1-p)]
HWE_3(p) = [p*p*p, 3*p*p*(1-p), 3*p*(1-p)*(1-p), (1-p)*(1-p)*(1-p)]
HWE_4(p) = [p*p*p*p, 4*p*p*p*(1-p), 6*p*p*(1-p)*(1-p), 4*p*(1-p)*(1-p)*(1-p), (1-p)*(1-p)*(1-p)*(1-p)]
HWE(p) = [HWE_1(p), HWE_2(p), HWE_3(p), HWE_4(p)]
end

# ╔═╡ 5ab3ba10-9955-11eb-1a12-7522af278a81
HWEv = HWE(0.5)

# ╔═╡ 96ed7e30-998c-11eb-2a78-0b90cc4522c6
genos = [[1.,0.],[2.,1.,0.],[3.,2.,1.,0.],[4.,3.,2.,1.,0.]]

# ╔═╡ 64f75740-998a-11eb-2b32-3be7ffc3b2c6
begin
	linear_map(v,α,d)=d*α*v.*(2/v[1])
	exponential_map(v,α,d)=α*v.*(2/v[1]).^d
end

# ╔═╡ 24d80f30-998d-11eb-3ae4-1fb8cf59a7db
phenos = linear_map.(genos,1.,1.)

# ╔═╡ d6951824-f8e2-494a-91fa-4171aa017508
q = reduce(vcat,HWEv)

# ╔═╡ 9ea6ebb0-998d-11eb-26ba-eb12b920bf01
begin
	pHWE25 = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="D=0.25")
	c1 = 1
	for (i,t) in (enumerate(linear_map.(genos,1.,0.25)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c1], ms=7.5)
		c1+= 1
		end
	end
	xlabel!("\$k\$")
	ylabel!("Trait mean")
end

# ╔═╡ a091a4b0-998d-11eb-1d55-3f0288c81feb
begin
	pHWE50 = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="D=0.50")
	c2 = 1
	for (i,t) in (enumerate(linear_map.(genos,1.,0.50)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c2], ms=7.5)
		c2+= 1
		end
	end
	xlabel!("\$k\$")
	ylabel!("Trait mean")
end

# ╔═╡ a1e044c0-998d-11eb-0053-cb74ee4562da
begin
	pHWE75 = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="D=0.75")
	c3 = 1
	for (i,t) in (enumerate(linear_map.(genos,1.,0.75)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c3], ms=7.5)
		c3+= 1
		end
	end
	xlabel!("\$k\$")
	ylabel!("Trait mean")
end

# ╔═╡ 308c314e-998b-11eb-0e06-6765aa9acd7d
begin
	pHWE100 = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="D=1")
	c4 = 1
	for (i,t) in (enumerate(linear_map.(genos,1.,1.)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c4], ms=7.5)
		c4+= 1
		end
	end
	xlabel!("\$k\$")
	ylabel!("Trait mean")
end

# ╔═╡ 5050db7e-2b37-4d9c-ae07-4836409eefe0
begin
	pHWE0e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="D=0")
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

# ╔═╡ 556426ae-e1f2-4987-beec-6a42490a5ce3
begin
	pHWE25e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="D=0.25")
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

# ╔═╡ 1467a6ea-ea7e-47cb-af3a-bf64824750e9
begin
	pHWE50e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="D=0.50")
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

# ╔═╡ 16fb1f40-a966-4920-952b-c5d366b7108f
begin
	pHWE75e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="D=0.75")
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

# ╔═╡ fa082dbf-dbb9-4dda-b34b-0b4ed56c5082
begin
	pHWE100e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400), title="D=1")
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

# ╔═╡ 81f637a0-c81f-4733-b3a5-f43482dff943
plot(pHWE25,pHWE50,pHWE75,pHWE100)

# ╔═╡ 1999d7bc-844c-4c83-88b0-a9f080d591ba
plot(pHWE0e,pHWE25e,pHWE50e,pHWE75e,pHWE100e)

# ╔═╡ 21c69f40-9876-11eb-0f98-cf5c73410b72
md""" ###### Effect of genotype to phenotype map"""

# ╔═╡ 9f0dcd88-a765-49fd-9b2c-73210a61618a
#add domninance, recessive

# ╔═╡ 06fdb4ff-0ae0-4cab-872f-248c57636bd0
t = randagent_p(0.5, 1., 100, [0., 0., 0., 4.],1)[1]

# ╔═╡ eb3861f9-82c8-42f5-96f0-8a5835dd2fea
t.loci

# ╔═╡ abe2c5fa-66c3-4596-8597-1e0f631ac9fa
trait_dom(t)

# ╔═╡ 525af8b0-c1ab-46b5-bb4d-608fb5001dc3
trait_rec(t)

# ╔═╡ 671e7783-258f-4ac8-854c-c74abb4e7c65
trait_exp(a::Agent) = (sum(a))*(ploidy(a)^(-a.d))

# ╔═╡ 7d5350c9-8e87-4488-9832-27de32f3fb97
begin
td2 = MixedPloidyDeme(agents = randagent_p(0.5, 1., 1, [0., 1., 0., 0.], 200), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 1., Vs = 1.)
td4 = MixedPloidyDeme(agents = randagent_p(0.5, 1., 1, [0., 0., 0., 1.], 200), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 1., Vs = 1.)
end

# ╔═╡ 268eb9a9-57e2-4a42-9b3f-f5866fbbfb2f
malthusian_fitness(td2,t,trait_add)

# ╔═╡ b5de9dfd-b875-4189-880e-d3b49a14d299
malthusian_fitness(td2,t,trait_dom)

# ╔═╡ bb6b43f2-d6f3-4322-9906-4e75698527a9
malthusian_fitness(td2,t,trait_rec)

# ╔═╡ 6a68dd7b-65d4-4f58-aff2-371d79f17064
begin
add_p2 = evolving_selectiondeme(td2,trait_add,500)
add_p4 = evolving_selectiondeme(td4,trait_add,500)
dom_p2 = evolving_selectiondeme(td2,trait_dom,500)
dom_p4 = evolving_selectiondeme(td4,trait_dom,500)
rec_p2 = evolving_selectiondeme(td2,trait_rec,500)
rec_p4 = evolving_selectiondeme(td4,trait_rec,500)
end

# ╔═╡ 13005ec3-a28c-4ac6-8cee-3d6c4a5def84
begin
plot(Normal(mean(dom_p2.fta[1]), std(dom_p2.fta[1])), color=:black, label="domp2")
plot!(Normal(mean(dom_p4.fta[1]), std(dom_p4.fta[1])), color=:blue, label="domp4")
plot!(Normal(mean(rec_p2.fta[1]), std(rec_p2.fta[1])), color=:green, label="rec2")
plot!(Normal(mean(rec_p4.fta[1]), std(rec_p4.fta[1])), color=:red, label="rec4")
plot!(Normal(mean(add_p2.fta[1]), std(add_p2.fta[1])), color=:orange, label="add2")
plot!(Normal(mean(add_p4.fta[1]), std(add_p4.fta[1])), color=:pink, label="add4")
xlabel!("Phenotype before selection")
end

# ╔═╡ 440b6029-d2f8-4f36-866e-9a79f68266bf
begin
	Hₒd_p2 = map(mean, dom_p2.het)
	Hₒd_p4 = map(mean, dom_p4.het)
	Hₒr_p2 = map(mean, rec_p2.het)
	Hₒr_p4 = map(mean, rec_p4.het)
	Hₒa_p2 = map(mean, add_p2.het)
	Hₒa_p4 = map(mean, add_p4.het)

	plot(Hₒd_p2, grid=false, color=:black, label="\$H_op1(t)\$", title="LOH")
	plot!(Hₒd_p4, grid=false, color=:blue, label="\$H_op2(t)\$")
	plot!(Hₒr_p2, grid=false, color=:red, label="\$H_op4(t)\$")
	plot!(Hₒr_p4, grid=false, color=:red, label="\$H_op4(t)\$")
	plot!(Hₒa_p2, grid=false, color=:red, label="\$H_op4(t)\$")
	plot!(Hₒa_p4, grid=false, color=:red, label="\$H_op4(t)\$")
	#plot!(0:sim_p2.ngen+1, 
		#t->expected_heterozygosity(0.5*(1-0.5), t, 50),
		#linestyle=:dash, color=:black, 
		#label = "\$(1-1/N)^t H_o(0)\$")
	plot!(0:neutral_p1.ngen+1, 
		t->eh_1(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:black, 
		label = "\$(1-1/N)^t H_o(0)\$")
	plot!(0:neutral_p2.ngen+1, 
		t->eh_2(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:blue, 
		label = "\$(1-1/2N)^t H_o(0)\$")
	plot!(0:neutral_p4.ngen+1, 
		t->eh_4(0.5*(1-0.5), t, 100),
		linestyle=:dash, color=:red, 
		label = "\$(1-1/4N)^t H_o(0)\$")
	xlabel!("\$t\$")
	ylabel!("\$H(t)\$")
end

# ╔═╡ 6030f17d-1ce9-42fe-9c90-44a4101f7312
begin
	traitmean_p2r = map(mean, rec_p2.tm)
	traitmean_p4r = map(mean, rec_p4.tm)
	traitmean_p2d = map(mean, dom_p2.tm)
	traitmean_p4d = map(mean, dom_p4.tm)
	traitmean_p2a = map(mean, add_p2.tm)
	traitmean_p4a = map(mean, add_p4.tm)
	plot(traitmean_p2r, grid=false, color=:black, label=false)
	plot!(traitmean_p4r, grid=false, color=:blue, label=false)
	plot!(traitmean_p2d, grid=false, color=:black, label=false)
	plot!(traitmean_p4d, grid=false, color=:blue, label=false)
	plot!(traitmean_p2a, grid=false, color=:black, label=false)
	plot!(traitmean_p4a, grid=false, color=:blue, label=false)
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 79b23d8f-18b3-416f-8a44-0d4baf3d9cec
begin
	popsize_p2r = map(mean, rec_p2.pop)
	popsize_p4r = map(mean, rec_p4.pop)
	popsize_p2d = map(mean, dom_p2.pop)
	popsize_p4d = map(mean, dom_p4.pop)
	popsize_p2a = map(mean, add_p2.pop)
	popsize_p4a = map(mean, add_p4.pop)
	plot(popsize_p2r, grid=false, color=:black, label=false)
	plot!(popsize_p4r, grid=false, color=:blue, label=false)
	plot!(popsize_p2d, grid=false, color=:blue, label=false)
	plot!(popsize_p4d, grid=false, color=:blue, label=false)
	plot!(popsize_p2a, grid=false, color=:blue, label=false)
	plot!(popsize_p4a, grid=false, color=:blue, label=false)
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ e40a4e8e-d06f-423a-bd7d-d9758cbfec48
begin
lin_p2 = evolving_selectiondeme(td2,trait_add,500)
lin_p4 = evolving_selectiondeme(td4,trait_add,500)
exp_p2 = evolving_selectiondeme(td2,trait_exp,500)
exp_p4 = evolving_selectiondeme(td4,trait_exp,500)
end

# ╔═╡ 6a699fb2-48fc-4218-aab4-df4b9ff24149
begin
#plot(Normal(mean(add_p2.fta[1]), std(add_p2.fta[1])), color=:orange, label="add2")
#plot!(Normal(mean(add_p4.fta[1]), std(add_p4.fta[1])), color=:pink, label="add4")
plot(Normal(mean(exp_p2.fta[1]), std(exp_p2.fta[1])), color=:blue, label="exp2")
plot!(Normal(mean(exp_p4.fta[1]), std(exp_p4.fta[1])), color=:red, label="exp4")
end

# ╔═╡ 9515122f-382b-4ffe-962a-be1d7e40526b
begin
plot(Normal(mean(lin_p2.fta[1]), std(lin_p2.fta[1])), color=:green, label="lin2")
plot!(Normal(mean(lin_p4.fta[1]), std(lin_p4.fta[1])), color=:pink, label="lin4")
end

# ╔═╡ 65b6bbe9-da25-4e77-a21d-196dd4d60df9
std(exp_p2.fta[1]), std(exp_p4.fta[1])

# ╔═╡ 6764af4b-03f1-4522-a29e-c707547774c2


# ╔═╡ cb9a459b-e6c1-4aaf-98c4-075ad8e79cc9
begin
td002 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 1., 0., 0.], 200, d=0.0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td022 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 1., 0., 0.], 200, d=0.2), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td042 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 1., 0., 0.], 200, d=0.4), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td062 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 1., 0., 0.], 200, d=0.6), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td082 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 1., 0., 0.], 200, d=0.8), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td102 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 1., 0., 0.], 200, d=1.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td004 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 0., 0., 1.], 200, d=0.), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td024 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 0., 0., 1.], 200, d=0.2), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td044 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 0., 0., 1.], 200, d=0.4), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td064 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 0., 0., 1.], 200, d=0.6), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td084 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 0., 0., 1.], 200, d=0.8), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
td104 = MixedPloidyDeme(agents = randagent_p(0.1, 1., 100, [0., 0., 0., 1.], 200, d=1.0), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 10., Vs = 1.)
end

# ╔═╡ 0bfb7600-2649-49e9-a5f3-bef66dbf3659
begin
exp_p002 = evolving_selectiondeme(td002,trait_exp,500)
exp_p102 = evolving_selectiondeme(td102,trait_exp,500)
exp_p004 = evolving_selectiondeme(td004,trait_exp,500)
exp_p024 = evolving_selectiondeme(td024,trait_exp,500)
exp_p044 = evolving_selectiondeme(td044,trait_exp,500)
exp_p064 = evolving_selectiondeme(td064,trait_exp,500)
exp_p084 = evolving_selectiondeme(td084,trait_exp,500)
exp_p104 = evolving_selectiondeme(td104,trait_exp,500)
end

# ╔═╡ 6b178069-09d4-42e8-85bc-8c36b81eb762
begin
plot(Normal(mean(exp_p002.fta[1]), std(exp_p002.fta[1])), color=:blue, label="exp002")
plot!(Normal(mean(exp_p102.fta[1]), std(exp_p102.fta[1])), color=:blue, label="exp102")
plot!(Normal(mean(exp_p004.fta[1]), std(exp_p004.fta[1])), color=:red, label="exp004")
plot!(Normal(mean(exp_p024.fta[1]), std(exp_p024.fta[1])), color=:red, label="exp024")
plot!(Normal(mean(exp_p044.fta[1]), std(exp_p044.fta[1])), color=:red, label="exp044")
plot!(Normal(mean(exp_p064.fta[1]), std(exp_p064.fta[1])), color=:red, label="exp064")
plot!(Normal(mean(exp_p084.fta[1]), std(exp_p084.fta[1])), color=:red, label="exp084")
plot!(Normal(mean(exp_p104.fta[1]), std(exp_p104.fta[1])), color=:red, label="exp104")

end

# ╔═╡ aa08d942-7a5c-40c5-ab08-bb3d9193f632
traitmean_p104 = map(mean, exp_p104.tm)

# ╔═╡ c1b4d088-5bc7-44f2-bc8b-f7dc34b5a61f
begin
	p104 = plot(traitmean_p104, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright, title="Tetraploid")
	for (i,t) in enumerate(exp_p104.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	#hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	#ylims!(18,22)
end

# ╔═╡ b1d74108-092a-448e-9154-7cee90312ae9
md""" #### Environmental stochasticity """

# ╔═╡ 7747be96-6a81-42f6-b9ff-c69ef44b4648
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

# ╔═╡ 1572237c-8151-4a73-b89c-696841a0b620
env = 0.75 .* sin.(1/50 .* (1:500)) .+ 20

# ╔═╡ 8b9268c8-e570-45c0-9fae-dc5c3dbd113e
plot(env)

# ╔═╡ 3b67e761-38ef-45ab-8ad3-76db4681bdd3
d_p2env = MixedPloidyDeme(agents = randagent_p(0.5, 0.1, 200, [0., 1., 0., 0.],1), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.95 0.05 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], θ = 20., Vs = 0.1)

# ╔═╡ 922118f6-2492-438a-b5c7-e0824b9ed0dd
stabselenv_p2 = evolving_selectiondemeenv(d_p2env, trait_exp, env, 500)

# ╔═╡ 694dce49-d920-4470-8767-44993a973cf9
begin
traitmean_p2env = map(mean, stabselenv_p2.tm)
popsize_p2env = map(mean, stabselenv_p2.pop)
end

# ╔═╡ 4f94e771-4b7d-4c5e-b71c-983fe829171b
begin
	p2env = plot(traitmean_p2env, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright, title="Diploid")
	for (i,t) in enumerate(stabselenv_p2.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(18,22)
end

# ╔═╡ 7b466623-5b18-4027-8829-424671044dd8
begin
	plot(stabselenv_p2.pop, grid=false, color=:black, label=false)
	plot!(stabselenv_p2.p2, grid=false, color=:red, label=false)
	plot!(stabselenv_p2.p4, grid=false, color=:blue, label=false)
	xlabel!("\$t\$")
	ylabel!("Population size")
end

# ╔═╡ Cell order:
# ╟─b79eb4e0-8696-11eb-0949-555c0c3c411f
# ╠═8ea08d1e-8bff-11eb-3fbb-f74cf53a2eed
# ╠═ded34c80-8699-11eb-33ce-214cb4f59699
# ╟─09e48960-9875-11eb-35fe-5122ae8dbe5f
# ╠═217c5170-8e52-11eb-175c-c9b8476ae81b
# ╠═3f70dcc0-8bfd-11eb-3df6-6709900e71ce
# ╠═39b5e750-8dbd-11eb-3f64-c12fa3fe01d5
# ╟─3106b1d0-9875-11eb-36ac-0d2b98a73fd4
# ╠═a8d89200-8e52-11eb-1277-b1857f7ea0fe
# ╠═66a25c90-8bff-11eb-2a53-09b98728f17c
# ╠═621f41a0-9885-11eb-2bf8-a530b5cc7e86
# ╠═a6806a20-9927-11eb-3fa6-3769de70a1f5
# ╠═3ec21690-992c-11eb-3b9f-4143af21edc5
# ╠═415089a2-992c-11eb-24c9-ff5fa01c0273
# ╠═42f31610-992c-11eb-03e1-5b5f8ca6a2fc
# ╠═1b1541de-8dbf-11eb-041c-815c116efeeb
# ╟─3fc09382-986b-11eb-29d6-57eeb0eddb77
# ╠═47dc60c0-9876-11eb-3d70-cd03509cbad8
# ╠═d6c03b40-9939-11eb-244b-3d747272f440
# ╠═fc96ac10-9884-11eb-1735-d751fd99c805
# ╠═66257e10-9879-11eb-2367-d706689f326a
# ╠═0501879e-9928-11eb-01c1-23fdad736664
# ╠═70a06c90-992a-11eb-1a39-ff907051a107
# ╠═7209e7a0-992a-11eb-0690-159af85e9c8d
# ╠═63ad1c50-987f-11eb-1512-df6fc858aaca
# ╠═44614280-9920-11eb-302b-6f959c8c2b94
# ╠═809cf8b2-9926-11eb-2b5e-d97d23d44531
# ╠═651f5aa2-9927-11eb-392f-c7944c0ad0aa
# ╠═134a9ebe-9921-11eb-3d7e-cdd5d3ff9ee1
# ╠═bae32f30-9926-11eb-3f23-cf639e204a0e
# ╠═c1e57720-9926-11eb-1d2d-058c7c88e98d
# ╠═158facc0-9921-11eb-25ef-a759b775ecdd
# ╠═cb25c830-9926-11eb-01f3-b3ddd767318b
# ╠═32ea0626-9019-4f27-8844-70c4f696549e
# ╟─db982b20-9923-11eb-285b-b9d57489b052
# ╠═c9b2b8d0-9879-11eb-294a-952184f8bf9c
# ╠═9cdac49e-9884-11eb-2be7-832739ca5d66
# ╟─fe5aa0f0-9939-11eb-1416-91041276c925
# ╠═5ca18050-9955-11eb-2703-47174bc221a8
# ╠═5ab3ba10-9955-11eb-1a12-7522af278a81
# ╠═96ed7e30-998c-11eb-2a78-0b90cc4522c6
# ╠═64f75740-998a-11eb-2b32-3be7ffc3b2c6
# ╠═24d80f30-998d-11eb-3ae4-1fb8cf59a7db
# ╠═d6951824-f8e2-494a-91fa-4171aa017508
# ╠═9ea6ebb0-998d-11eb-26ba-eb12b920bf01
# ╠═a091a4b0-998d-11eb-1d55-3f0288c81feb
# ╠═a1e044c0-998d-11eb-0053-cb74ee4562da
# ╠═308c314e-998b-11eb-0e06-6765aa9acd7d
# ╠═5050db7e-2b37-4d9c-ae07-4836409eefe0
# ╠═556426ae-e1f2-4987-beec-6a42490a5ce3
# ╠═1467a6ea-ea7e-47cb-af3a-bf64824750e9
# ╠═16fb1f40-a966-4920-952b-c5d366b7108f
# ╠═fa082dbf-dbb9-4dda-b34b-0b4ed56c5082
# ╠═81f637a0-c81f-4733-b3a5-f43482dff943
# ╠═1999d7bc-844c-4c83-88b0-a9f080d591ba
# ╟─21c69f40-9876-11eb-0f98-cf5c73410b72
# ╠═9f0dcd88-a765-49fd-9b2c-73210a61618a
# ╠═91826552-6b17-491f-ad1b-77e06dfcdeb6
# ╠═06fdb4ff-0ae0-4cab-872f-248c57636bd0
# ╠═eb3861f9-82c8-42f5-96f0-8a5835dd2fea
# ╠═abe2c5fa-66c3-4596-8597-1e0f631ac9fa
# ╠═525af8b0-c1ab-46b5-bb4d-608fb5001dc3
# ╠═671e7783-258f-4ac8-854c-c74abb4e7c65
# ╠═7d5350c9-8e87-4488-9832-27de32f3fb97
# ╠═268eb9a9-57e2-4a42-9b3f-f5866fbbfb2f
# ╠═b5de9dfd-b875-4189-880e-d3b49a14d299
# ╠═bb6b43f2-d6f3-4322-9906-4e75698527a9
# ╠═6a68dd7b-65d4-4f58-aff2-371d79f17064
# ╠═13005ec3-a28c-4ac6-8cee-3d6c4a5def84
# ╟─440b6029-d2f8-4f36-866e-9a79f68266bf
# ╠═6030f17d-1ce9-42fe-9c90-44a4101f7312
# ╠═79b23d8f-18b3-416f-8a44-0d4baf3d9cec
# ╠═e40a4e8e-d06f-423a-bd7d-d9758cbfec48
# ╠═6a699fb2-48fc-4218-aab4-df4b9ff24149
# ╠═9515122f-382b-4ffe-962a-be1d7e40526b
# ╠═65b6bbe9-da25-4e77-a21d-196dd4d60df9
# ╠═6764af4b-03f1-4522-a29e-c707547774c2
# ╠═cb9a459b-e6c1-4aaf-98c4-075ad8e79cc9
# ╠═0bfb7600-2649-49e9-a5f3-bef66dbf3659
# ╠═6b178069-09d4-42e8-85bc-8c36b81eb762
# ╠═aa08d942-7a5c-40c5-ab08-bb3d9193f632
# ╠═c1b4d088-5bc7-44f2-bc8b-f7dc34b5a61f
# ╠═b1d74108-092a-448e-9154-7cee90312ae9
# ╠═7747be96-6a81-42f6-b9ff-c69ef44b4648
# ╠═1572237c-8151-4a73-b89c-696841a0b620
# ╠═8b9268c8-e570-45c0-9fae-dc5c3dbd113e
# ╠═f4182f59-e44c-4f1d-838b-f4cd9ae8f22d
# ╠═3b67e761-38ef-45ab-8ad3-76db4681bdd3
# ╠═922118f6-2492-438a-b5c7-e0824b9ed0dd
# ╠═694dce49-d920-4470-8767-44993a973cf9
# ╠═4f94e771-4b7d-4c5e-b71c-983fe829171b
# ╠═7b466623-5b18-4027-8829-424671044dd8
