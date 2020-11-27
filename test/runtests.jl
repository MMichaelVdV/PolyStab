using Test
using PolyStab: SimpleDeme, randagent, random_mating, allelefreqs

@testset "haploid drift" begin
    d = SimpleDeme(agents=randagent(0.2, 0.1, 100, 100))
    d = foldl((d,i)->random_mating(d), 1:1000, init=d)
    @test sort(unique(allelefreqs(d))) == [0., 1.] 
end
