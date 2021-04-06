using Test
#using PolyStab: SimpleDeme, randagent, random_mating, allelefreqs
using PolyStab: MixedPloidyDeme, randagent, randagent_p, random_mating, allelefreqs_p, mate_p, trait, mean, var



#@testset "haploid drift" begin
#    d = SimpleDeme(agents=randagent(0.2, 0.1, 100, 100))
#    d = foldl((d,i)->random_mating(d), 1:1000, init=d)
#    @test sort(unique(allelefreqs(d))) == [0., 1.] 
#end

@testset "Polyploid drift" begin
    dp = MixedPloidyDeme(agents = randagent_p(0.5, 1., 200, [0., 1., 0., 0.], 100, d=2.), OV = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [0. 0. 0. 0. ; 1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.] )
    dp = foldl((dp,i)->random_mating(dp), 1:2500, init=dp)
	@test sort(unique(allelefreqs_p(dp))) == [0., 1.]
end

@testset "Segregation" begin
    n = 1000
    k = 2
    α = 1/√(k*n) 
    #α = 1.
    a = randagent(0.5, α, n, k)
    b = randagent(0.5, α, n, k)
    o = map(x->trait(mate_p(a,b)), 1:2500) 
    @info "" mean(o) var(o)
    @test mean(o) ≈ 0.5*(trait(a) + trait(b)) atol=0.1
    @test var(o) ≈ n * (α^2 / 4) atol=0.1 
end