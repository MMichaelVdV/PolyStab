#### PolyStab: Individual-based modeling to study polyploid establishment from an eco-evo perspective

This repository contains an individual-based model and accompanying Pluto notebooks with toy simulations 
to study polyploid establishment from an integrated eco-evolutionary persepective. 

##### Overview of the Pluto notebooks:

- unreduced_gametes.jl: Simulations with mixed-ploidy populations in a single deme
- AM.jl: Single deme simulations with assortative mating and/or selfing
- scaling_phenotype.jl: More in depth look at genotype-phenotype map
- gen_var.jl: Simulations of genetic variance in a single deme and different genotype-phenotype maps
- Env_stoch.jl: Changing environment through time in a single deme
- FToNS.jl: Testing some different model parameters in a single deme
- fixedpop.jl: Selection/Truncsel with constant population size
- truncation_selection.jl: Some experiments with truncation selection in a single deme
- island.jl: Implementation of island migration model
- island_estab_sims.jl: Simulations for polyploid establishment with island migration
- model_with_migration.jl: Original PnB simulations with haploid population
- linhab_gifs.jl: PnB type simulations with mixed ploidy populations in a linear habitat
- habitat2D.jl: Preliminary simulations in 2D habitat

Old code (might still contain some notes):
- matrix.jl
- migration_gifs.jl
- mixed_ploidy.jl
- single_deme_selection.jl
- single_deme.jl

