# How to compute the scale factors for the energy responce correction in FastSim ECAL

## Generation of Particle Gun samples (without correction)

Insert the requested energy in `test/SingleElectron_fast_template_cfi_GEN_SIM.py`
and `test/SingleElectron_full_template_cfi_GEN_SIM.py` and generate some events.

Optional: run `test/prepareFastSim.py` to run many jobs with different energies.

## Save the imporant information in a tree

Extract the important information (genEta, genEnergy, response) from the generated samples
and storem them in a TTree, using `plugin/ConfFile_cfg.py` (here, you have to modify the input files as well)

In the end, you need two TTrees, one for fastsim, and one for fullsim

## Calculate the scales

Use the `test/unbinnedScaling.cc` program (compile with the given Makefile).

