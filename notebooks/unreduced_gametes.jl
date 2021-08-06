### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ d4f68ac6-3953-4040-8eff-d80fdc865ee6
begin
using GLM
using CSV
end

# ╔═╡ cc1dcf70-7c6c-11eb-3664-b7e70e49723a
using Parameters, Random, Distributions, Plots, StatsBase, PlutoUI, ColorSchemes

# ╔═╡ e419ec50-87fa-11eb-230a-8f6396f0e459
using Turing, MCMCChains, StatsPlots, DataFrames

# ╔═╡ 0144f6d0-87fb-11eb-0fa7-7149bfb316db
using StatsFuns: logistic

# ╔═╡ 47183362-97e7-11eb-2573-99a97f304e99
using PolyStab: Agent, randagent_p, MixedPloidyDeme, trait, evolving_ugdeme, evolving_selectiondeme, heterozygosities_p, allelefreqs_p, evolving_neutraldeme, recombination, random_mating, allelefreqs_p, heterozygosities_p, AbstractDeme, ploidy

# ╔═╡ bd2ab750-7c60-11eb-3cc5-6d57dd448dfd
md"""##### Minority Cytotype Exclusion in Local Plant Populations (Levin, 1995)"""

# ╔═╡ 98180cf0-7c6b-11eb-0da9-c55b3cf2a066


# ╔═╡ 2dadbbf0-7c5f-11eb-3d09-1fcf2b82ead9
md"""##### Establishment of a tetraploid cytotype in a diploid population: Effect of relative fitness of the cytotypes (Felber, 1991)"""

# ╔═╡ 472102d2-84a8-11eb-2d6d-eb101b7a707b


# ╔═╡ 7df71be0-7c70-11eb-3b92-91457ecf38fd
md"""### Case 1 (Felber, 1991)"""

# ╔═╡ c8734230-7c65-11eb-30bd-4314f44ee3ba
md"""  
- Fertility of diploids and tetraploids is the same
- Viability of diploids and tetraploids is the same
- Triploids are not viable
- (The system is at equilibrium)

The only parameter in this model is u: unreduced gamete formation for diploids.

Starting from a population of only diploid individuals, according to this model, tetraploids will exclude diploids if unreduced gamete formation for is higher than 0.1716.

"""

# ╔═╡ 575c04e0-97e7-11eb-0750-eb108125e90a
d_p = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],200), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.83 0.17 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.], K=200)

# ╔═╡ e0789de0-7c6e-11eb-0900-3387c904ea0b
sim_popvar = evolving_selectiondeme(d_p,150)

# ╔═╡ e5916552-7c6e-11eb-1868-a7b856e0ad1f
begin
	pf1_var = sim_popvar.p2
	pf2_var = sim_popvar.p4
	plot(pf1_var, grid=false, color=:blue, label="diploid", title="u=0.17")
	plot!(pf2_var, grid=false, color=:red, label="tetraploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ 4cde5ba0-7c79-11eb-23e1-d1ac938c0e6b
function random_search(t)
	OV = OffspringViability([[1.,0.,0.1,0.],[0.,1.,0.,0.],[0.1,0.,0.,0.],[0.,0.,0.,0.]])
	ploidy = []
	param = []
	for i in 1:t
		u = rand(Uniform(0., 0.5))
		UG = UnreducedGamete([[1-u,u,0.,0.],[0.,0.,0.,0.],[0.0,1.,0.0,0.]])
		sim_ploidyvar = evolving_deme_ploidyvar(d_p,100,UG,OV)
		if sim_ploidyvar.p2[end] >= sim_ploidyvar.p4[end]
			push!(ploidy,2)
		else
			push!(ploidy,4)
		end
		push!(param, u)
	end
	ploidy, param
end
		

# ╔═╡ 1f6a4e72-7d48-11eb-0933-eb25e7a5b568
md"""I used grid search on an interval for from 0 to 0.5 with 200 steps. For each value of u I simulated 50 generations and checked the most frequent ploidy level at the final generation"""

# ╔═╡ e315b920-7c7d-11eb-0f34-5782accf7286
function grid_search(t)
	ploidy = []
	param = []
	pop_size = []
	p2 = []
	p4 = []
	for u in range(0, stop=0.5, length=t)
		for rep in 1:10
		UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1. 0. 0.]
		d_p = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],45), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UG, K=50, Vs=1.2)
		sim_ploidyvar = evolving_selectiondeme(d_p,50)
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

# ╔═╡ 39c92860-7c7e-11eb-1efe-13bf99e99d1e
stats_2 = grid_search(100)

# ╔═╡ 9e59d2ef-7ced-497b-b938-52d20ac0a186
md""" ###### Fitting logit model"""

# ╔═╡ c88e13d0-7c7e-11eb-0119-731f7ddc65ce
dp_2 = [(stats_2[2][i],stats_2[1][i],stats_2[3][i]) for i in 1:1000]

# ╔═╡ f3d8c63c-2ae3-4fa1-91e1-48b08d71599b
begin
df = DataFrame([stats_2[1] stats_2[2]])
rename!(df,:x1 => :Ploidy)
rename!(df,:x2 => :u)
end

# ╔═╡ 5aa68d72-3533-4451-b2a0-c144cd287321
#begin
#	df = DataFrame(CSV.File("./model_1.csv"))
#	first(df,5)
#end

# ╔═╡ 876c4cba-eedc-4477-8e4c-3e94a1360d48
Y = Int.((df[1]./2).-1) 

# ╔═╡ b7874cc0-ca2a-4e3a-a926-8add52f71b54
df[1] = Y

# ╔═╡ 37231930-7c7f-11eb-2524-9115fc85d926
begin
p3 = scatter(dp_2, label=false, title="UG for 2n only")
vline!([0.17], label="u=0.17",linewidth=5)
xlabel!("u")
ylabel!("Ploidy")
end

# ╔═╡ ff0cb820-7ff9-11eb-00e7-cf90fe869537
begin
p7 = plot(stats_2[2], stats_2[3], grid=false, color=:white, label="Pop size after t generations")
scatter!(stats_2[2], stats_2[4], grid=false, color=:green, label="Diploids")
scatter!(stats_2[2], stats_2[5], grid=false, color=:red, label="Tetraploids")
vline!([0.17], label="u=0.17",linewidth=5)
xlabel!("\$u\$")
ylabel!("Number of individuals")
end

# ╔═╡ ec061ade-8801-11eb-2449-f3788b4ffa9b
md"""
#Bayesian logistic regression 

begin
M = [stats_2[1] stats_2[2]];
data = DataFrame(M, ["Ploidy","u"]);
using CSV
CSV.write("./model_1.csv", data)

Yz = Int.((stats_2[1]./2).-1) 
X = stats_2[2]
@model logreg(data, X) = begin
    n = 200
	m = 1
    σ ~ truncated(Cauchy(0, 1), 0., Inf)
	#intercept ~ Normal(0, σ)
    β ~ MvNormal(zeros(m), σ)
    Z = logistic.(X.*β) #intercept .+ 
    for i=1:n
        data[i] ~ Bernoulli(Z[i])
    end 
end
chain = sample(logreg(Y, X), NUTS(), 500) 	
plot(chain)
β = hcat(get(chain, :β).β...)
βmean = vec(mean(β, dims=1))

logistic.(X.*βmean)
end
"""

# ╔═╡ a478224e-7c70-11eb-133b-03e027bf81c8
md""" ### Case 2 (Felber, 1991)"""

# ╔═╡ 57d84bf0-7c6b-11eb-0e68-f99c717f46e4
md"""  
See case 1 but now the parameter in this model is u(=v): unreduced gamete formation for both diploids and tetraploids.

Starting from a population of only diploid individuals, according to this model, tetraploids will exclude diploids if unreduced gamete formation for is higher than 0.2.

"""

# ╔═╡ 0cdd11b0-97ed-11eb-1862-a1c728804ef5
d_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],175), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 0.80 0.20 0. 0. ; 0. 0. 0. 0. ; 0. 0.8 0. 0.2], K=200)

# ╔═╡ f89b70d0-7c70-11eb-2d9c-d9abf241dd17
sim_ploidyvar2 = evolving_selectiondeme(d_p2,500)

# ╔═╡ 28e30f50-7c71-11eb-072c-dffe12b4b437
begin
	pf2_p2 = sim_ploidyvar2.p2
	pf3_p2 = sim_ploidyvar2.p3
	pf4_p2 = sim_ploidyvar2.p4
	plot(pf2_p2, grid=false, color=:blue, label="diploid", title="u=v=0.1")
	#plot!(pf3_p2, grid=false, color=:green, label="triploid")
	plot!(pf4_p2, grid=false, color=:red, label="tetraploid")
	xlabel!("\$t\$")
	ylabel!("Number of individuals")
end

# ╔═╡ ed0024e0-7cd5-11eb-1548-f70db6fdda06
function grid_search_2(t)
	ploidy = []
	param = []
	pop_size = []
	p2 = []
	p4 = []
	for u in range(0, stop=0.5, length=t)
		for rep in 1:10
		UG = [0. 0. 0. 0. ; 1-u u 0. 0. ; 0. 0. 0. 0. ; 0. 1-u 0. u ]
		d_p2 = MixedPloidyDeme(agents = randagent_p(0.5, 0.5, 50, [0., 1., 0., 0.],175), OV = [1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = UG, K=200)
		sim_ploidyvar = evolving_selectiondeme(d_p2,50)
		if sim_ploidyvar.p2[end] >= sim_ploidyvar.p4[end]
			push!(ploidy,2)
		else
			push!(ploidy,4)
		end
		pop = (sim_ploidyvar.p2[end] + sim_ploidyvar.p4[end])
		push!(pop_size,pop)
		push!(param, u)
		push!(p2, sim_ploidyvar.p2[end])
		push!(p4, sim_ploidyvar.p4[end])
		end
	end
	ploidy, param, pop_size, p2, p4
end

# ╔═╡ 52fc8ae0-7cd6-11eb-1731-f3c71b08b7a9
stats_3 = grid_search_2(200)

# ╔═╡ 5d85451e-f189-44eb-9a24-463b865b0e97
md""" Remark: Probability of tetraploid establishment is here caclulated as ratio of simulations where a majority of the individuals of the population are tetraploids after 50 generation to the total number all simulations, and not conditional on a certain population size. This is why P estab starts to drop from a certain value of `u`, since populations will no longer be sustainable due to cytotype load and the total number of established populations will drop (both 2n and 4n). (Could also plot the ratio of 4n to 2n conditional on a certain population size higher than 0)"""

# ╔═╡ 99dcebd2-7cd6-11eb-3b96-85bc3a0940d8
dp_3 = [(stats_3[2][i],stats_3[1][i],stats_3[3][i]) for i in 1:2000]

# ╔═╡ c0370950-7cd6-11eb-1113-9b2c9104da16
begin
p4 = scatter(dp_3, label=false, title="UG for both 2n and 4n")
vline!([0.2], label="u=0.2",linewidth=5)
xlabel!("u")
ylabel!("Ploidy")
end

# ╔═╡ 55373330-7ff8-11eb-351c-9df644942e4d
begin
p8 = plot(stats_3[2], stats_3[3], grid=false, color=:white, label="Pop size after t generations")
scatter!(stats_3[2], stats_3[4], grid=false, color=:green, label="Diploids")
scatter!(stats_3[2], stats_3[5], grid=false, color=:red, label="Tetraploids")
vline!([0.17], label="u=0.17",linewidth=5)
xlabel!("\$u\$")
ylabel!("Number of individuals")
end

# ╔═╡ 1abbabe0-87fb-11eb-271e-53e02f1ec547
md""" ### Cytogenic load"""

# ╔═╡ 86be28f0-84ad-11eb-2ccf-6f60ad079d75
begin
#credits to arzwa
#This is the simple Felber model, expressed in terms of marginal/relative fitnesses
w = [0., 1., 0., 1.]
# wij = proportion of offspring of ploidy j for a given individual of ploidy i
# marginal/relative fitness of ploidy i (proportion of offspring in the next generation coming from a parent of ploidy i, regardless of offspring ploidy level) is
# wi* = ∑ⱼwij. The mean fitness is then ∑ᵢwi* pᵢ. The proportion of offspring
# of ploidy level j in the next generation is (∑ᵢwij pᵢ) / w̄. This seems like a
# mathematically quite elegant formulation of the model, and one that may
# generalize relatively easily.
w22(w, u, p) = w[2] * (1-u)^2 * p
w24(w, u, p) = w[4] * (u^2 * p + u*(1-u) * (1-p))
w44(w, u, p) = w[4] * ((1-u)^2 * (1-p))+ u*(1-u) * p #u*(1-u) * p: where does this come from ?
w̄(w, u, p) = w22(w, u, p) * p  + w44(w, u, p) * (1-p) + w24(w, u, p) * p     
function w̄_(w, u, p) 
    w2 = w22(w, u, p) * p  
    w4 = w44(w, u, p) * (1-p) + w24(w, u, p) * p 
    return w2, w4
end
function evolve(w, u, p, n)
    map(1:n) do x
        w2, w4 = w̄_(w, u, p)
        p = w2 / (w2 + w4)
        p, w2 + w4
    end
end
# In this model (no triploids, tetraploids and diploids equally fit a priori,
# shared unreduced gamete formation rate, no higher polyploids) tetraploids
# take over irrespective of the initial frequencies as soon as u > 0.2. Note
# though that values of u > 0.2 lead to a very serious load in the population,
# as at least a fraction (1-0.8^2) of inviable offspring will be generated
# each generation!

p = plot(grid=false, legend=false, xlabel="\$t\$", ylabel="\$\\bar{w}\$")
cols = cgrad(:magma, 0.05:0.01:0.45)
for u=0.05:0.01:0.45
    mw = last.(evolve(w, u, 1., 30))
    plot!(p, mw, color=cols[u], linewidth=2, alpha=0.5)
end
plot(p, size=(750,750))
#savefig("felber-wbar.pdf")
# Interestingly, mean fitness decreases in this system. Of course it is well
# known that in systems with frequency dependence mean fitness does not
# necessarily increase, but somehow I never checked this for the Felber model.
end

# ╔═╡ 7f32ce70-85ab-11eb-1c30-1ffd857a2545
evolve([0., 1., 0., 1.], 0.2, 1., 20)

# ╔═╡ f641bd10-85aa-11eb-03fc-0de54af93d35
begin
pf = plot(grid=false, legend=true, xlabel="\$t\$", ylabel="\$\\bar{w}\$")
#cols = cgrad(:magma, 0.05:0.01:0.45)
#for u=0.05:0.01:0.45
	u = 0.25
	s = (evolve(w, u, 1., 30))
    mw = last.(s)
	pfs = first.(s)
	#f = 
    plot!(pf, mw, color=cols[u], linewidth=2, alpha=0.5, label="mean fitness")
	plot!(pf, pfs, color="blue", linewidth=2, alpha=0.5, label="diploid freq")
	plot!(pf, 1 .-pfs, color="red", linewidth=2, alpha=0.5, label="tetraploid freq")
#end
plot(pf, size=(750,750))
#savefig("felber-fubar.pdf")
# Interestingly, mean fitness decreases in this system. Of course it is well
# known that in systems with frequency dependence mean fitness does not
# necessarily increase, but somehow I never checked this for the Felber model.
end

# ╔═╡ c164fb61-7701-452b-8b8d-d6af89f7c091
fm = @formula(Ploidy ~ u)

# ╔═╡ 9cc7a396-fff2-4c50-85bf-ef75c574d975
logit = glm(fm, df, Binomial(), LogitLink())

# ╔═╡ 8a935580-85ad-11eb-1c01-fbfe6528c0e2
(evolve(w, u, 1., 30))

# ╔═╡ 495309e0-85ac-11eb-3231-b3dce7ae61c5
first.(evolve([0., 1., 0., 1.], 0.5, 1., 30))

# ╔═╡ aab09870-84b1-11eb-11c6-ffa1224114d8
begin
#Here I try to plot the cytogenic load for different levels of u to check if there is indeed a maximimum for a mixed population
cytoload_2n(u) = 2*(1-u)*u
cytoload_4n(u) = 2*(1-u)*u+u^2
#cytoload_total(u,f) = (2*(1-u)*u)*(1-f) + (2*(1-u)*u+u^2)*(f)
#as total frequency of unviable gametes
cytoload_total(u,f) = (1-u)*f*(1-u)*(1-f)+(1-u)*f*u*(1-f)+(1-u)*f*u*f+(1-u)*(1-f)*u*(1-f)+(u*(1-f))^2+u*f*u*(1-f)
	
f=[i for i in 0.0:0.01:1.]

pl = plot(grid=false, legend=false, xlabel="\$f(4N)\$", ylabel="\$cytoload\$")
#cols = cgrad(:magma, 0.05:0.01:0.45)
for u=0.05:0.01:0.45
	load = cytoload_total.(u,f)
    plot!(f, load, color=cols[u], linewidth=2, alpha=0.5)
end
plot(pl, size=(750,750))
	
end

# ╔═╡ 302363ee-7c6c-11eb-0d7b-effeadcdeb60
md"""#### Functions"""

# ╔═╡ b54e9160-98c0-11eb-0f5a-f1216dc4f615
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

# ╔═╡ a352f65e-98c3-11eb-2048-d3e731947c4b
begin
tick1 = stabprob(stats_2[1])
grid1 = plot([0.005:0.005:0.5...],tick1,label=false)
vline!([0.17],label="u=0.17",linewidth=2,style=:dash)
hline!([0.50],label=false,linewidth=2,style=:dash)
xlabel!("u")
ylabel!("P estab")

it(x) = 1/(1+exp(-(73.1085*x-10.4124)))
vline!([10.4124/73.1085], linewidth=2,style=:dash, label="u_crit")
plot!(df[2],it.(df[2]), colour =:black, label=false)
end

# ╔═╡ f5064130-98c1-11eb-1d9c-bbe9ff236b23
begin
tick2 = stabprob(stats_3[1])
grid2 = plot([0.0025:0.0025:0.5...],tick2, label=false)
vline!([0.20],label="u=0.20",linewidth=2,style=:dash)
hline!([0.50],label=false,linewidth=2,style=:dash)
xlabel!("u")
ylabel!("P estab")
end

# ╔═╡ d75aee00-7c84-11eb-2d5a-b3ca3c4bd625
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
	up2, up4
end

# ╔═╡ 41392210-7c85-11eb-19bf-675f7a929063
binz = bin(dp_2)

# ╔═╡ b7e83d6e-7c84-11eb-0e72-7b4c25184221
begin
	p1 = histogram(binz[1], bins=50, label=false, title="Diploid > Tetraploid")
	vline!([0.17], label="u=0.17",linewidth=5)
	xlabel!("u")
	ylabel!("Count")
end

# ╔═╡ 616389e0-7c85-11eb-13c9-fbf0587d3ccc
begin
	p2 = histogram(binz[2], bins=50, label=false, title="Tetraploid > Diploid")
	vline!([0.17], label="u=0.17",linewidth=5)
	xlabel!("u")
	ylabel!("Count")
end

# ╔═╡ 7f182630-7cc6-11eb-071d-0d22322e64e2
plot(p1,p2,grid1,p7, legend=false)

# ╔═╡ fb49b790-7cd6-11eb-3043-3563217174e9
binz2 = bin(dp_3)

# ╔═╡ 08b7378e-7cd7-11eb-0120-edc0d8f9dbb8
begin
	p5 = histogram(binz2[1], bins=50, label=false, title="Diploid > Tetraploid")
	vline!([0.2], label="u=0.2",linewidth=5)
	xlabel!("u=v")
	ylabel!("Count")
end

# ╔═╡ 38ef6cc0-7cd7-11eb-26e6-1767c49d50de
begin
	p6 = histogram(binz2[2], bins=50, label=false, title="Tetraploid > Diploid")
	vline!([0.2], label="u=0.2",linewidth=5)
	xlabel!("u=v")
	ylabel!("Count")
end

# ╔═╡ c5025790-7cd7-11eb-24f8-51d85dc4b968
plot(p5,p6,grid2,p8, legend=false)

# ╔═╡ Cell order:
# ╟─bd2ab750-7c60-11eb-3cc5-6d57dd448dfd
# ╠═98180cf0-7c6b-11eb-0da9-c55b3cf2a066
# ╟─2dadbbf0-7c5f-11eb-3d09-1fcf2b82ead9
# ╠═472102d2-84a8-11eb-2d6d-eb101b7a707b
# ╟─7df71be0-7c70-11eb-3b92-91457ecf38fd
# ╟─c8734230-7c65-11eb-30bd-4314f44ee3ba
# ╠═575c04e0-97e7-11eb-0750-eb108125e90a
# ╠═e0789de0-7c6e-11eb-0900-3387c904ea0b
# ╠═e5916552-7c6e-11eb-1868-a7b856e0ad1f
# ╟─4cde5ba0-7c79-11eb-23e1-d1ac938c0e6b
# ╟─1f6a4e72-7d48-11eb-0933-eb25e7a5b568
# ╠═e315b920-7c7d-11eb-0f34-5782accf7286
# ╠═39c92860-7c7e-11eb-1efe-13bf99e99d1e
# ╟─9e59d2ef-7ced-497b-b938-52d20ac0a186
# ╠═c88e13d0-7c7e-11eb-0119-731f7ddc65ce
# ╠═f3d8c63c-2ae3-4fa1-91e1-48b08d71599b
# ╠═d4f68ac6-3953-4040-8eff-d80fdc865ee6
# ╠═5aa68d72-3533-4451-b2a0-c144cd287321
# ╠═876c4cba-eedc-4477-8e4c-3e94a1360d48
# ╠═b7874cc0-ca2a-4e3a-a926-8add52f71b54
# ╠═c164fb61-7701-452b-8b8d-d6af89f7c091
# ╠═9cc7a396-fff2-4c50-85bf-ef75c574d975
# ╠═a352f65e-98c3-11eb-2048-d3e731947c4b
# ╠═37231930-7c7f-11eb-2524-9115fc85d926
# ╠═41392210-7c85-11eb-19bf-675f7a929063
# ╠═b7e83d6e-7c84-11eb-0e72-7b4c25184221
# ╠═616389e0-7c85-11eb-13c9-fbf0587d3ccc
# ╠═ff0cb820-7ff9-11eb-00e7-cf90fe869537
# ╠═7f182630-7cc6-11eb-071d-0d22322e64e2
# ╠═ec061ade-8801-11eb-2449-f3788b4ffa9b
# ╟─a478224e-7c70-11eb-133b-03e027bf81c8
# ╟─57d84bf0-7c6b-11eb-0e68-f99c717f46e4
# ╠═0cdd11b0-97ed-11eb-1862-a1c728804ef5
# ╠═f89b70d0-7c70-11eb-2d9c-d9abf241dd17
# ╟─28e30f50-7c71-11eb-072c-dffe12b4b437
# ╠═ed0024e0-7cd5-11eb-1548-f70db6fdda06
# ╠═52fc8ae0-7cd6-11eb-1731-f3c71b08b7a9
# ╠═f5064130-98c1-11eb-1d9c-bbe9ff236b23
# ╟─5d85451e-f189-44eb-9a24-463b865b0e97
# ╠═99dcebd2-7cd6-11eb-3b96-85bc3a0940d8
# ╠═c0370950-7cd6-11eb-1113-9b2c9104da16
# ╠═55373330-7ff8-11eb-351c-9df644942e4d
# ╠═fb49b790-7cd6-11eb-3043-3563217174e9
# ╟─08b7378e-7cd7-11eb-0120-edc0d8f9dbb8
# ╟─38ef6cc0-7cd7-11eb-26e6-1767c49d50de
# ╠═c5025790-7cd7-11eb-24f8-51d85dc4b968
# ╟─1abbabe0-87fb-11eb-271e-53e02f1ec547
# ╠═86be28f0-84ad-11eb-2ccf-6f60ad079d75
# ╟─7f32ce70-85ab-11eb-1c30-1ffd857a2545
# ╠═f641bd10-85aa-11eb-03fc-0de54af93d35
# ╟─8a935580-85ad-11eb-1c01-fbfe6528c0e2
# ╟─495309e0-85ac-11eb-3231-b3dce7ae61c5
# ╟─aab09870-84b1-11eb-11c6-ffa1224114d8
# ╟─302363ee-7c6c-11eb-0d7b-effeadcdeb60
# ╠═cc1dcf70-7c6c-11eb-3664-b7e70e49723a
# ╠═e419ec50-87fa-11eb-230a-8f6396f0e459
# ╠═0144f6d0-87fb-11eb-0fa7-7149bfb316db
# ╠═47183362-97e7-11eb-2573-99a97f304e99
# ╠═b54e9160-98c0-11eb-0f5a-f1216dc4f615
# ╠═d75aee00-7c84-11eb-2d5a-b3ca3c4bd625
