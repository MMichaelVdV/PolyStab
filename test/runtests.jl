using Test
#using PolyStab: SimpleDeme, randagent, random_mating, allelefreqs
using PolyStab: SimpleDeme, MixedPloidyDeme, randagent, randagent_p, random_matingh, random_mating, mating_PnB 
using PolyStab: allelefreqs_p, mate_p, trait, trait_mean, mean, var

@testset "Drift" begin
    d = SimpleDeme(agents=randagent_p(0.5, 1., 200, [1., 0., 0., 0.], 100))
    d = foldl((d,i)->random_matingh(d), 1:2500, init=d)
    @test sort(unique(allelefreqs_p(d))) == [0., 1.] 

    dp = MixedPloidyDeme(agents=randagent_p(0.5, 0.125, 200, [0., 1., 0., 0.], 100), OV=[1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG=[0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.] )
    dp = foldl((dp,i)->random_mating(dp), 1:2500, init=dp)
	@test sort(unique(allelefreqs_p(dp))) == [0., 1.]
end

@testset "Selection" begin
    ds = MixedPloidyDeme(agents=randagent_p(0.5, 0.125, 200, [0., 1., 0., 0.], 100), OV=[1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG=[0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.] )
    ds = foldl((dp,i)->mating_PnB(dp), 1:2500, init=ds)
    @info "" trait_mean(ds) length(ds)
    @test sort(unique(allelefreqs_p(ds))) == [0., 1.]
    @test trait_mean(ds) ≈ ds.θ atol=ds.θ/10
    @test length(ds) ≈ ds.K atol=ds.K/5
    #have to set atol range broad enough because of stochasticity or better to take mean of last n generations?
end

@testset "Segregation" begin
    n = 50
    k = 2
    #α = 1/√(k*n) 
    α = .5
    a = randagent(0.5, α, n, k)
    b = randagent(0.5, α, n, k)
    o = map(x->trait(mate_p(a,b)), 1:2500) 
    @info "" mean(o) var(o) 0.5*(trait(a)+trait(b)) n*(α^2 / 4)
    @test mean(o) ≈ 0.5*(trait(a) + trait(b)) atol=0.1
    @test var(o) ≈ 0.5*n * (1/k)*(α^2 / 4) atol=0.1
end