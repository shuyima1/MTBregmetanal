function [flux] = BiomassComponentFluxes(model)

biomasscomponentRxns = {'PROT';'RNA';'DNA';'SM_MOL';'PE';'TAG';'PIMS';'LAM';'MAPC';'P-L-GLX'};
soln = optimizeCbModel(model);

flux = [soln.f;soln.x];