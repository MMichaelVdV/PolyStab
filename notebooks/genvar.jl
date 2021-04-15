### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ 8ea08d1e-8bff-11eb-3fbb-f74cf53a2eed
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes, DataFrames, StatsPlots

# ╔═╡ ded34c80-8699-11eb-33ce-214cb4f59699
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, evolving_neutraldeme, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy, trait_mean, randagent, evolving_haploiddeme, mate_p, evolving_selectiondemeh

# ╔═╡ b79eb4e0-8696-11eb-0949-555c0c3c411f
md"""### Genetic variance in a mixed ploidy population"""

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
neutral_p1 = evolving_haploiddeme(d_p1, 500)
neutral_p2 = evolving_neutraldeme(d_p2, 500)
neutral_p4 = evolving_neutraldeme(d_p4, 500)
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

O_mean2 = mean(o2)
E_mean2 = 0.5*(trait(a2) + trait(b2)) 
O_var2 = var(o2)
E_var2 = n * (α^2 / 4)
	
a4 = randagent(0.5, α, n, 4)
b4 = randagent(0.5, α, n, 4)
o4 = map(x->trait(mate_p(a4,b4)), 1:5000) 

O_mean4 = mean(o4)
E_mean4 = 0.5*(trait(a4) + trait(b4)) 
O_var4 = var(o4)
E_var4 = 
	
plot(Normal(O_mean2, sqrt(O_var2)),color=:blue, label="Observed 2N")
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
	p1n = plot(traitneutral_p1, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright)
	for (i,t) in enumerate(neutral_p1.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(18,22)
end

# ╔═╡ 415089a2-992c-11eb-24c9-ff5fa01c0273
begin
	p2n = plot(traitneutral_p2, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright)
	for (i,t) in enumerate(neutral_p2.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(18,22)
end

# ╔═╡ 42f31610-992c-11eb-03e1-5b5f8ca6a2fc
begin
	p4n = plot(traitneutral_p4, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright)
	for (i,t) in enumerate(neutral_p4.fta)
	scatter!([i for x in 1:10],t,label=false,colour="black",ma=0.35,ms=2.5)
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
	hline!([d_p1.θ],label="Optimal phenotype",colour="black",linestyle=:dash)
	ylims!(18,22)
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
stabsel_p1 = evolving_selectiondemeh(d_p1,500) #need to fix, switches to diploid
stabsel_p2 = evolving_selectiondeme(d_p2,500)
stabsel_p4 = evolving_selectiondeme(d_p4,500)
end

# ╔═╡ d6c03b40-9939-11eb-244b-3d747272f440
begin
plot(Normal(mean(stabsel_p1.fta[1]), std(stabsel_p1.fta[1])), color=:black, label="Haploid")
plot!(Normal(mean(stabsel_p2.fta[1]), std(stabsel_p2.fta[1])), color=:blue, label="Diploid")
plot!(Normal(mean(stabsel_p4.fta[1]), std(stabsel_p4.fta[1])), color=:red, label="Tetraploid")
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
	p1 = plot(traitmean_p1, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright)
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
	p2 = plot(traitmean_p2, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright)
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
	p4 = plot(traitmean_p4, grid=false, color=:red, label=false,linewidth=3,legend=:bottomright)
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
	plot(Normal(mean(stabsel_p1.fta[1]), std(stabsel_p1.fta[1])), label="")
	plot!(Normal(mean(stabsel_p1.fta[end]), std(stabsel_p1.fta[end])), label="")
end

# ╔═╡ 809cf8b2-9926-11eb-2b5e-d97d23d44531
histogram(stabsel_p1.fta[1], bins = 25, fillalpha = 0.4)

# ╔═╡ 651f5aa2-9927-11eb-392f-c7944c0ad0aa
histogram(stabsel_p1.fta[end], bins = 25, fillalpha = 0.4)

# ╔═╡ 134a9ebe-9921-11eb-3d7e-cdd5d3ff9ee1
begin
	plot(Normal(mean(stabsel_p2.fta[1]), std(stabsel_p2.fta[1])), label="")
	plot!(Normal(mean(stabsel_p2.fta[end]), std(stabsel_p2.fta[end])), label="")
end

# ╔═╡ bae32f30-9926-11eb-3f23-cf639e204a0e
histogram(stabsel_p2.fta[1], bins = 25, fillalpha = 0.4)

# ╔═╡ c1e57720-9926-11eb-1d2d-058c7c88e98d
histogram(stabsel_p2.fta[end], bins = 25, fillalpha = 0.4)

# ╔═╡ 158facc0-9921-11eb-25ef-a759b775ecdd
begin
	plot(Normal(mean(stabsel_p4.fta[1]), std(stabsel_p4.fta[1])), label="")
	plot!(Normal(mean(stabsel_p4.fta[end]), std(stabsel_p4.fta[end])), label="")
end

# ╔═╡ cb25c830-9926-11eb-01f3-b3ddd767318b
histogram(stabsel_p4.fta[1], bins = 25, fillalpha = 0.4)

# ╔═╡ d0d5eda0-9926-11eb-30ac-9998afb1a028
histogram(stabsel_p4.fta[end], bins = 25, fillalpha = 0.4)

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
	pHWE25 = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400))
	c1 = 1
	for (i,t) in (enumerate(linear_map.(genos,1.,0.25)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c1], ms=7.5)
		c1+= 1
		end
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ a091a4b0-998d-11eb-1d55-3f0288c81feb
begin
	pHWE50 = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400))
	c2 = 1
	for (i,t) in (enumerate(linear_map.(genos,1.,0.50)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c2], ms=7.5)
		c2+= 1
		end
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ a1e044c0-998d-11eb-0053-cb74ee4562da
begin
	pHWE75 = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400))
	c3 = 1
	for (i,t) in (enumerate(linear_map.(genos,1.,0.75)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c3], ms=7.5)
		c3+= 1
		end
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 308c314e-998b-11eb-0e06-6765aa9acd7d
begin
	pHWE100 = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400))
	c4 = 1
	for (i,t) in (enumerate(linear_map.(genos,1.,1.)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c4], ms=7.5)
		c4+= 1
		end
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 5050db7e-2b37-4d9c-ae07-4836409eefe0
begin
	pHWE0e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400))
	c5 = 1
	for (i,t) in (enumerate(exponential_map.(genos,1.,0.)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c5], ms=7.5)
		c5+= 1
		end
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 556426ae-e1f2-4987-beec-6a42490a5ce3
begin
	pHWE25e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400))
	c6 = 1
	for (i,t) in (enumerate(exponential_map.(genos,1.,0.25)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c6], ms=7.5)
		c6+= 1
		end
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 1467a6ea-ea7e-47cb-af3a-bf64824750e9
begin
	pHWE50e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400))
	c7 = 1
	for (i,t) in (enumerate(exponential_map.(genos,1.,0.50)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c7], ms=7.5)
		c7+= 1
		end
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 16fb1f40-a966-4920-952b-c5d366b7108f
begin
	pHWE75e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400))
	c8 = 1
	for (i,t) in (enumerate(exponential_map.(genos,1.,0.75)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c8], ms=7.5)
		c8+= 1
		end
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ fa082dbf-dbb9-4dda-b34b-0b4ed56c5082
begin
	pHWE100e = plot(grid=false, color=:black,label=false,linewidth=3,legend=:false, size=(400,400))
	c10 = 1
	for (i,t) in (enumerate(exponential_map.(genos,1.,1.)))
		for x in 1:length(t)
	scatter!(([i],t[x]),label=false,colour="black",ma=q[c10], ms=7.5)
		c10+= 1
		end
	end
	xlabel!("\$t\$")
	ylabel!("Trait mean")
end

# ╔═╡ 81f637a0-c81f-4733-b3a5-f43482dff943
plot(pHWE25,pHWE50,pHWE75,pHWE100)

# ╔═╡ 1999d7bc-844c-4c83-88b0-a9f080d591ba
plot(pHWE0e,pHWE25e,pHWE50e,pHWE75e,pHWE100e)

# ╔═╡ 21c69f40-9876-11eb-0f98-cf5c73410b72
md""" ###### Effect of genotype to phenotype map"""

# ╔═╡ 9f0dcd88-a765-49fd-9b2c-73210a61618a
#add domninance

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
# ╟─fc96ac10-9884-11eb-1735-d751fd99c805
# ╠═66257e10-9879-11eb-2367-d706689f326a
# ╠═0501879e-9928-11eb-01c1-23fdad736664
# ╠═70a06c90-992a-11eb-1a39-ff907051a107
# ╠═7209e7a0-992a-11eb-0690-159af85e9c8d
# ╟─63ad1c50-987f-11eb-1512-df6fc858aaca
# ╠═44614280-9920-11eb-302b-6f959c8c2b94
# ╠═809cf8b2-9926-11eb-2b5e-d97d23d44531
# ╠═651f5aa2-9927-11eb-392f-c7944c0ad0aa
# ╠═134a9ebe-9921-11eb-3d7e-cdd5d3ff9ee1
# ╠═bae32f30-9926-11eb-3f23-cf639e204a0e
# ╠═c1e57720-9926-11eb-1d2d-058c7c88e98d
# ╠═158facc0-9921-11eb-25ef-a759b775ecdd
# ╠═cb25c830-9926-11eb-01f3-b3ddd767318b
# ╠═d0d5eda0-9926-11eb-30ac-9998afb1a028
# ╟─db982b20-9923-11eb-285b-b9d57489b052
# ╟─c9b2b8d0-9879-11eb-294a-952184f8bf9c
# ╟─9cdac49e-9884-11eb-2be7-832739ca5d66
# ╟─fe5aa0f0-9939-11eb-1416-91041276c925
# ╟─5ca18050-9955-11eb-2703-47174bc221a8
# ╟─5ab3ba10-9955-11eb-1a12-7522af278a81
# ╟─96ed7e30-998c-11eb-2a78-0b90cc4522c6
# ╠═64f75740-998a-11eb-2b32-3be7ffc3b2c6
# ╠═24d80f30-998d-11eb-3ae4-1fb8cf59a7db
# ╟─d6951824-f8e2-494a-91fa-4171aa017508
# ╟─9ea6ebb0-998d-11eb-26ba-eb12b920bf01
# ╟─a091a4b0-998d-11eb-1d55-3f0288c81feb
# ╟─a1e044c0-998d-11eb-0053-cb74ee4562da
# ╟─308c314e-998b-11eb-0e06-6765aa9acd7d
# ╠═5050db7e-2b37-4d9c-ae07-4836409eefe0
# ╟─556426ae-e1f2-4987-beec-6a42490a5ce3
# ╟─1467a6ea-ea7e-47cb-af3a-bf64824750e9
# ╠═16fb1f40-a966-4920-952b-c5d366b7108f
# ╟─fa082dbf-dbb9-4dda-b34b-0b4ed56c5082
# ╠═81f637a0-c81f-4733-b3a5-f43482dff943
# ╠═1999d7bc-844c-4c83-88b0-a9f080d591ba
# ╟─21c69f40-9876-11eb-0f98-cf5c73410b72
# ╠═9f0dcd88-a765-49fd-9b2c-73210a61618a
