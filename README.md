# WaspSim: a MATLAB simulation model to evaluate biocontrol feasibility for invasive insects
The WaspSim model is a MATLAB package developed as part of a project aimed at evaluating whether an invasive insect would be a good candidate for a classical biological control programme. The model could help agencies avoid the time and expense of full economic evaluations for projects that are technically doomed.

For details of the model and an application see:

Cacho, O.J. and Hester, S.M. (2022). Modelling biocontrol of invasive insects: an application to European Wasp (Vespula germanica) in Australia. Ecological Modeling. 467, 109939.           https://doi.org/10.1016/j.ecolmodel.2022.109939.

To use the model download all the files and subfolders from the Matlab folder. There is no preparation required to run the existing application. All the data required to run the application is available in Matlab (.mat) files and will load as the code is executed.

To get started, open Sim_cluster_script.m in the Matlab environment. This script allows you to generate the results reported in Cacho and Hester (2022). Before running the code, consider that it may take several days to complete the simulations, depending on the computer used. You may wish to change the experimental design and simulate only a few cases for a single spatial cluster to get a feel for the model.

The Data folder contains four results files  nsw_20_10_1.mat, vic_20_10_1.mat, nsw_50_01_1.mat and vic_50_01_1.mat which are used by Plot_figs.m to produce phase diagrams and other useful plots without having to first run Sim_cluster_script.m. 

