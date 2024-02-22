# HyperonProduction

HyperonProduction is a LArSoft analyser module targeting sigma zero hyperon
selections using MicroBooNE MC and data files.

As part of a larger LArSoft build: this is meant to be placed within a version
of `ubana` tagged as `v08_00_00_43`. The module is configured through the
`job/hyperonConfig.fcl` and can be invoked through using
`job/run_HyperonProduction.fcl` with the `lar` command.

The analyser module operates on `reco2` files that have been processed using
the standard/near-standard MicroBooNE production workflow.
