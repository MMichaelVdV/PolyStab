using Test
using PolyStab: SimpleDeme, randagent, random_mating, allelefreqs

@testset "haploid drift" begin
    d = SimpleDeme(agents=randagent(0.2, 0.1, 100, 100))
    d = foldl((d,i)->random_mating(d), 1:1000, init=d)
    @test sort(unique(allelefreqs(d))) == [0., 1.] 
end

@testset "polyploid drift" begin
    dp = MixedPloidyDeme(agents = vcat(randagent_p(0.5, 1., 200, [2], 100, d=2.),randagent_p(0.5, 0.1, 250, [4], 0, d=4.)), OV = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.], UG = [1. 0. 0. 0. ; 0. 0. 0. 0. ; 0. 0. 0. 0.] )
    dp = foldl((dp,i)->random_mating_mixedp(dp), 1:1000, init=dp)
	sort(unique(allelefreqs_p(dp))) == [0., 1.]