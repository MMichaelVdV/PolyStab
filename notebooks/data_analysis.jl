### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 246b9d1f-3f4b-46a3-9ac7-3cd33d728922
begin
using Pkg
using DataFrames
using CSV
using Plots
using GLM
using StatsBase
using Lathe
using MLBase
using ClassImbalance
using ROCAnalysis
end

# ╔═╡ eda9be70-04fc-4b14-9257-ccace4fccabd
using PyCall

# ╔═╡ 7091d208-4840-484f-b97b-bb849fab2646
using Lathe.preprocess: TrainTestSplit

# ╔═╡ 6fb96da4-3aa8-4e5c-8756-00e64e1fbec2
begin
using Turing, Distributions
using RDatasets
using MCMCChains, StatsPlots
using StatsFuns: logistic
using MLDataUtils: shuffleobs, stratifiedobs, rescale!
using Random
Random.seed!(0)
Turing.setprogress!(false)
end

# ╔═╡ 5b2317eb-0e42-49d2-92f2-df3a5d1b38f6
md"""### Logistic regression: the classical approach"""

# ╔═╡ 98d665f0-e8be-11eb-0559-5b90358cc7c6
begin
	df = DataFrame(CSV.File("./model_1.csv"))
	first(df,5)
end

# ╔═╡ c6d1550f-a12d-4e4c-ac7f-ce0f1d4ef68a
size(df)

# ╔═╡ 2dba49ba-0147-49ee-b684-00cc4f723b6f
describe(df)

# ╔═╡ 0f3e85dd-8986-434f-b742-1804ea37ea05
names(df)

# ╔═╡ 9f723833-3d9a-4250-a48f-27e4fd624bd2
countmap(df.Ploidy)

# ╔═╡ 5b7fe268-35d5-4452-a5cb-3d4f871befda
Y = Int.((df[1]./2).-1) 

# ╔═╡ 5d2c988d-9abe-41d9-aad4-dd4014aa8c23
df[1] = Y

# ╔═╡ ef400ff9-1664-4e81-bd63-96b361b5b4d3
df

# ╔═╡ 16033ca0-80ce-4642-904f-02aaefb15bba
df[1]

# ╔═╡ de9833de-5b9e-407a-bc11-e35e83cb242e
fm = @formula(Ploidy ~ u)

# ╔═╡ 09ae38b0-ccf7-4b1d-9721-031964a7cf59
#train, test = TrainTestSplit(df,1.)

# ╔═╡ 9087d7b8-3020-4d0e-84c6-75effee6af6a
logit = glm(fm, df, Binomial(), LogitLink())

# ╔═╡ ec9632e7-c2e0-4ecb-ba5e-82a2e2fbb20a
it(x) = 1/(1+exp(-(28.7974*x-4.35138)))

# ╔═╡ c3ebcce0-8f03-42b1-a0ef-e51ba07331ff
plot(df[2],it.(df[2]), colour =:black, label=false)

# ╔═╡ 5d674532-e227-4ab1-8043-637f12df44ee


# ╔═╡ 5bc87b8a-a93d-42b7-abef-c67c67762598
md"""### Logistic regression: the Bayesian approach"""

# ╔═╡ 84cafaf6-80e3-4cec-9f99-1d635b1a004f
data = df

# ╔═╡ 99a21c74-fdb8-407f-929b-73fe0de51d67
begin
	# Convert "Default" and "Student" to numeric values.
	df[!,:Ploidy] = [r.Ploidy == 4 ? 1.0 : 0.0 for r in eachrow(data)]
	
	# Delete the old columns which say "Yes" and "No".
	#select!(data, Not([:Default, :Student]))
	
	# Show the first six rows of our edited dataset.
	first(df, 6)
end

# ╔═╡ 6652f82c-2067-4414-9d34-81ba5125ea62
begin
	function split_data(df, target; at = 0.70)
	    shuffled = shuffleobs(df)
	    trainset, testset = stratifiedobs(row -> row[target], 
	                                      shuffled, p = at)
	end
	
	features = [:u]
	numerics = [:u]
	target = :Ploidy
	
	trainset, testset = split_data(data, target, at = 0.05)
	for feature in numerics
	  μ, σ = rescale!(trainset[!, feature], obsdim=1)
	  rescale!(testset[!, feature], μ, σ, obsdim=1)
	end
	
	# Turing requires data in matrix form, not dataframe
	train = Matrix(trainset[:, features])
	test = Matrix(testset[:, features])
	train_label = trainset[:, target]
	test_label = testset[:, target];
end

# ╔═╡ 1d32d773-edc7-46f3-9108-d25c0f8981c2


# ╔═╡ 0dde0297-3e9d-4626-8d5d-459e0205c0da
# Bayesian logistic regression (LR)
begin
@model logistic_regression(x, y, n, σ) = begin
    intercept ~ Normal(0, σ)

    u ~ Normal(0, σ)
    

    for i = 1:n
        v = logistic(intercept + u*x[i, 1])
        y[i] ~ Bernoulli(v)
    end
end
	
# Retrieve the number of observations.
n, _ = size(train)

# Sample using HMC.
chain = mapreduce(c -> sample(logistic_regression(train, train_label, n, 1), HMC(0.05, 10), 1500),
    chainscat,
    1:3
)

describe(chain)
end

# ╔═╡ Cell order:
# ╟─5b2317eb-0e42-49d2-92f2-df3a5d1b38f6
# ╠═246b9d1f-3f4b-46a3-9ac7-3cd33d728922
# ╠═eda9be70-04fc-4b14-9257-ccace4fccabd
# ╠═98d665f0-e8be-11eb-0559-5b90358cc7c6
# ╠═c6d1550f-a12d-4e4c-ac7f-ce0f1d4ef68a
# ╠═2dba49ba-0147-49ee-b684-00cc4f723b6f
# ╠═99a21c74-fdb8-407f-929b-73fe0de51d67
# ╠═0f3e85dd-8986-434f-b742-1804ea37ea05
# ╠═9f723833-3d9a-4250-a48f-27e4fd624bd2
# ╠═7091d208-4840-484f-b97b-bb849fab2646
# ╠═5b7fe268-35d5-4452-a5cb-3d4f871befda
# ╠═5d2c988d-9abe-41d9-aad4-dd4014aa8c23
# ╠═ef400ff9-1664-4e81-bd63-96b361b5b4d3
# ╠═16033ca0-80ce-4642-904f-02aaefb15bba
# ╠═de9833de-5b9e-407a-bc11-e35e83cb242e
# ╠═09ae38b0-ccf7-4b1d-9721-031964a7cf59
# ╠═9087d7b8-3020-4d0e-84c6-75effee6af6a
# ╠═ec9632e7-c2e0-4ecb-ba5e-82a2e2fbb20a
# ╠═c3ebcce0-8f03-42b1-a0ef-e51ba07331ff
# ╠═5d674532-e227-4ab1-8043-637f12df44ee
# ╟─5bc87b8a-a93d-42b7-abef-c67c67762598
# ╠═6fb96da4-3aa8-4e5c-8756-00e64e1fbec2
# ╠═84cafaf6-80e3-4cec-9f99-1d635b1a004f
# ╠═6652f82c-2067-4414-9d34-81ba5125ea62
# ╠═1d32d773-edc7-46f3-9108-d25c0f8981c2
# ╠═0dde0297-3e9d-4626-8d5d-459e0205c0da
