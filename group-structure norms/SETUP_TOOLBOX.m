%Set up the different paths


%Path to SDPT3
path_to_sdpt3 = % *** TO BE FILLED ***


addpath(pwd);
addpath([pwd,'/Auxiliary_functions/']);
addpath([pwd,'/SDPT3_optimization/']);
addpath([pwd,'/Mex/']);
addpath([pwd,'/Experiments/']);

addpath([path_to_sdpt3,'/']);
addpath([path_to_sdpt3,'/Solver']);
addpath([path_to_sdpt3,'/Solver/Mexfun']);
addpath([path_to_sdpt3,'/Linsysolver/spchol']);
addpath([path_to_sdpt3,'/HSDSolver']);
addpath([path_to_sdpt3,'/Examples']);
addpath([path_to_sdpt3,'/testdir']);


% Compile mex files
mex Mex/mexGetFringeNonzeroPatterns.cpp
mex Mex/mexAndLogic.cpp
mex Mex/mexMultiColon.cpp

mex Mex/mexEta_update_for_regression.cpp
mex Mex/mexEta_update_for_regression_squared_penalty.cpp
mex Mex/mexW_Regularized_update.cpp -l blas -l lapack