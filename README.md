# acousticSimWrapper
Wrapper for performing an acoustic simulation using k-wave

## Getting Started

 * clone project and get k-wave library
 * add k-wave folder to matlab path and run the example_handheldsimulation.m

for lightfluence simulations
 * get our diffusive lightfluence simulation (currently not open source) but any other package will do as well
 * simulate the initial pressure distribution p_0 (in most of the examples it's a 300 x 300 x #wavelength matrix)

## Prerequisites - run with

 * [k-wave](httP://www.ka-wave.org)
 * for diffusive light-fluence simulation: getLightSourcesHH.m, getFluenceModel.m, make_fvm_mesh.m, make_fvm_forward.m
 * for light simulation: absorption spectra of Hb, HbO2, fat, water as described in spectra/spectrafiles.txt (they can be downloaded from http://omlc.org/spectra/hemoglobin/index.html)

## creating own simulations

 * AcousticSim()

## running GPU examples

to run simulation codes more efficient, the GPU support of k-wave should be used.
 * download kspaceFirstOrder3D-CUDA.exe and put it into the kwavegpusims folder

