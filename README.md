# HyperonProduction

HyperonProduction is a LArSoft analyser module targeting sigma zero hyperon
selections using MicroBooNE MC and data files.

The analyser module operates on `reco2` files that have been processed using
the standard/near-standard MicroBooNE production workflow. The primary output
of this module is processed ROOT ntuple data containing information on particle
interactions within the MicroBooNE TPC. Two ntuples are created: one with
per-event information simulated (if available) and reconstructed information,
and another with auxiliary information stored per sub-run period.

Partial documentation and a description of the format of output files can be
viewed at the [HyperonProduction
Reference](https://playonverbs.github.io/sigmazerosearch/larsoft/) and
[HyperonProduction NTuple
Specification](https://playonverbs.github.io/sigmazerosearch/ntuples/) pages.

As part of a larger LArSoft build: this is meant to be placed within a version
of `ubana` tagged as `v08_00_00_85`. The module is configured through the
`job/hyperonConfig.fcl` and can be invoked through using
`job/run_HyperonProduction*.fcl` with the `lar` command.
