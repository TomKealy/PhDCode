function launcher(vec)

% This is the script that has been used for comparing the computational
% efficiency of ISTA, and FISTA.
% it has to be launcher with toyBench.m


addpath('../ksvd/');
addpath('../ksvd/');
addpath('../ksvd/ompbox/');
addpath('../ksvd/tools/');
addpath('..')

%--------------------------------------------------------------------------
%           Problem parameters
%--------------------------------------------------------------------------
learning.dim=100;        % problem dimension
learning.TailleDic=10;   % Dictionary size
learning.T=3;            % nb de fonctions de bases
learning.K=200;           % nb de signaux
learning.nbitermax=30;
learning.rsnr=3;
learning.verbose=1;
learning.type='gaussian'; % is changed in constant or smooth
learning.nbtest=5000;
learning.nbval=200;
learning.lambdadico=0;
comment='';

%--------------------------------------------------------------------------
%       Dictionary Learning Constraint Options
%--------------------------------------------------------------------------
tvbound=0.1;
lambdasparseunitnorm=0.05;



%--------------------------------------------------------------------------
%       Dictionary Learning Options
%--------------------------------------------------------------------------

options.nbitermax=50;                   % nb of iterations
options.threshold=1e-4;                 % threshold on relative variation
options.warmstart=0;
options.eta=1;                          % multiplicative decrease of stopping
% criterion
options.ALnbitermax=100;
options.ALmu=10;
options.ALmudico=1;
%--------------------------------------------------------------------------
%  Majorization method options (Iterative Thresholding algorithm)
%--------------------------------------------------------------------------
options.it.threshold=1e-6;
options.it.nbitermax=2000;

%--------------------------------------------------------------------------
% Fista / Augmented Lagrangian options
%--------------------------------------------------------------------------
% matrix approximation options
options.matapprox.nbitermax=1000;
options.matapprox.threshold=1e-8;  % avec tv faut prendre + petit
options.matapprox.thresholdmin=1e-8;
% regularizer options
options.matapprox.rightregularizer='\ell_1-\ell_1';
options.matapprox.rightregparam=0;

%-------------------------------------------------------------------------
%           dictionary options
%--------------------------------------------------------------------------





options.matapprox.leftregularizer='admmvect';
operators(1).type='tv';
operators(1).parameters.bound=1e-2;   % a diminuer si moins de bruit
operators(1).parameters.nbitertv=300;
operators(2).type='unit-norm';
optionsadmm.nbitermax=500;
optionsadmm.seuil=1e-5;
options.matapprox.leftregparam.optionsadmm=optionsadmm;
options.matapprox.leftregparam.operatorsadmm=operators;


for i=1:length(vec)
    switch vec(i)
        case 1
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]
            learning.type='gaussian';
            comment='unit-norm';
            options.matapprox.leftregularizer='proj-unitnorm2';
            toyBenchMSE(comment,learning,options);
        case 2
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
           learning.lambdasparsecodevec=[0.05]
            learning.type='smooth';
            comment='unit-norm';
            options.matapprox.leftregularizer='proj-unitnorm2';
            toyBenchMSE(comment,learning,options);
        case 3
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]
            learning.type='constant';
            comment='unit-norm';
            options.matapprox.leftregularizer='proj-unitnorm2';
            toyBenchMSE(comment,learning,options);

        case 4
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]
            learning.type='constant';
            comment=['tvunit-norm'];
            options.comment.tv=num2str(tvbound);

            options.matapprox.leftregularizer='admmvect';
            operators(1).type='tv';
            operators(1).parameters.bound=tvbound;   % a diminuer si moins de bruit
            operators(1).parameters.nbitertv=300;
            operators(1).parameters.threshold=1e-5;
            operators(1).parameters.warmstart=1;
            operators(2).type='unit-norm';
            optionsadmm.nbitermax=500;
            optionsadmm.seuil=1e-5;
            options.matapprox.leftregparam.optionsadmm=optionsadmm;
            options.matapprox.leftregparam.operatorsadmm=operators;
            toyBenchMSE(comment,learning,options);
        case 5
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]
            learning.type='smooth';
            comment='tvunit-norm';
            options.comment.tv=num2str(tvbound);
            options.matapprox.leftregularizer='admmvect';
            operators(1).type='tv';
            operators(1).parameters.bound=tvbound;   % a diminuer si moins de bruit
            operators(1).parameters.nbitertv=300;
            operators(1).parameters.threshold=1e-5;
            operators(1).parameters.warmstart=1;
            operators(2).type='unit-norm';
            optionsadmm.nbitermax=500;
            optionsadmm.seuil=1e-5;
            options.matapprox.leftregparam.optionsadmm=optionsadmm;
            options.matapprox.leftregparam.operatorsadmm=operators;
            toyBenchMSE(comment,learning,options);
        case 6
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]
            learning.type='constant';
            comment='tv_posunit-norm';
            options.comment.tv=num2str(tvbound);
            options.matapprox.leftregularizer='admmvect';
            operators(1).type='tv';
            operators(1).parameters.bound=tvbound;   % a diminuer si moins de bruit
            operators(1).parameters.nbitertv=300;
            operators(1).parameters.threshold=1e-5;
            operators(1).parameters.warmstart=1;
            operators(2).type='unit-norm';
            operators(3).type='positive';
            optionsadmm.nbitermax=500;
            optionsadmm.seuil=1e-5;
            options.matapprox.leftregparam.optionsadmm=optionsadmm;
            options.matapprox.leftregparam.operatorsadmm=operators;
            toyBenchMSE(comment,learning,options);
        case 7
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]
            learning.type='smooth';
            comment='tv_posunit-norm';
            options.comment.tv=num2str(tvbound);
            options.matapprox.leftregularizer='admmvect';
            operators(1).type='tv';
            operators(1).parameters.bound=tvbound;   % a diminuer si moins de bruit
            operators(1).parameters.nbitertv=300;
            operators(1).parameters.threshold=1e-5;
            operators(1).parameters.warmstart=1;
            operators(2).type='unit-norm';
            operators(3).type='positive';

            optionsadmm.nbitermax=500;
            optionsadmm.seuil=1e-5;
            options.matapprox.leftregparam.optionsadmm=optionsadmm;
            options.matapprox.leftregparam.operatorsadmm=operators;
            toyBenchMSE(comment,learning,options);
        case 8
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]; % regularization parameters
            learning.type='smooth';
            comment='tv_sparseunit-norm';
            options.comment.tv=num2str(tvbound);
            options.comment.sparse=num2str(lambdasparseunitnorm);
            options.matapprox.leftregularizer='admmvect';
            operators(1).type='tv';
            operators(1).parameters.bound=tvbound;   % a diminuer si moins de bruit
            operators(1).parameters.nbitertv=300;
            operators(1).parameters.threshold=1e-5;
            operators(1).parameters.warmstart=1;
            operators(2).type='sparse_unit-norm';
            operators(2).parameters.threshold=lambdasparseunitnorm;
            optionsadmm.nbitermax=500;
            optionsadmm.seuil=1e-5;
            options.matapprox.leftregparam.optionsadmm=optionsadmm;
            options.matapprox.leftregparam.operatorsadmm=operators;
            toyBenchMSE(comment,learning,options);
        case 9
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]; % regularization parameters
            learning.type='constant';
            comment='tv_sparseunit-norm';
            options.comment.tv=num2str(tvbound);
            options.comment.sparse=num2str(lambdasparseunitnorm);
            options.matapprox.leftregularizer='admmvect';
            operators(1).type='tv';
            operators(1).parameters.bound=tvbound;   % a diminuer si moins de bruit
            operators(1).parameters.nbitertv=300;
            operators(1).parameters.threshold=1e-5;
            operators(1).parameters.warmstart=1;
            operators(2).type='sparse_unit-norm';
            operators(2).parameters.threshold=lambdasparseunitnorm;
            optionsadmm.nbitermax=500;
            optionsadmm.seuil=1e-5;
            options.matapprox.leftregparam.optionsadmm=optionsadmm;
            options.matapprox.leftregparam.operatorsadmm=operators;
            toyBenchMSE(comment,learning,options);
        case 10
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]; % regularization parameters
            learning.type='constant';
            comment='sparseunit-norm';
            lambdasparseunitnorm=0.1;
            options.comment.sparse=num2str(lambdasparseunitnorm);
            options.matapprox.leftregularizer='admmvect';
            operators(1).type='sparse_unit-norm';
            operators(1).parameters.threshold=lambdasparseunitnorm;
            optionsadmm.nbitermax=500;
            optionsadmm.seuil=1e-5;
            options.matapprox.leftregparam.optionsadmm=optionsadmm;
            options.matapprox.leftregparam.operatorsadmm=operators;
            toyBenchMSE(comment,learning,options);
        case 11
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]; % regularization parameters
            learning.type='constant';
            comment='posunit-norm';
            %options.comment.tv=num2str(tvbound);
            %options.comment.sparse=num2str(lambdasparseunitnorm);
            options.matapprox.leftregularizer='admmvect';
            operators(1).type='unit-norm';
            operators(2).type='positive';
            optionsadmm.nbitermax=500;
            optionsadmm.seuil=1e-5;
            options.matapprox.leftregparam.optionsadmm=optionsadmm;
            options.matapprox.leftregparam.operatorsadmm=operators;
            toyBenchMSE(comment,learning,options);
        case 12
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]; % regularization parameters
            learning.type='smooth';
            comment='sparseunit-norm';
            lambdasparseunitnorm=0.1;
            options.comment.sparse=num2str(lambdasparseunitnorm);
            options.matapprox.leftregularizer='admmvect';
            operators(1).type='sparse_unit-norm';
            operators(1).parameters.threshold=lambdasparseunitnorm;
            optionsadmm.nbitermax=500;
            optionsadmm.seuil=1e-5;
            options.matapprox.leftregparam.optionsadmm=optionsadmm;
            options.matapprox.leftregparam.operatorsadmm=operators;
            toyBenchMSE(comment,learning,options);
        case 13
            learning.lambdasparsecodevec=[0.01 0.05 0.1 0.5]; % regularization parameters
            learning.lambdasparsecodevec=[0.05]; % regularization parameters
            learning.type='constant';
            comment='tv_sparseunit-normpos';
            options.comment.tv=num2str(tvbound);
            options.comment.sparse=num2str(lambdasparseunitnorm);
            options.matapprox.leftregularizer='admmvect';
            operators(1).type='tv';
            operators(1).parameters.bound=tvbound;   % a diminuer si moins de bruit
            operators(1).parameters.nbitertv=300;
            operators(1).parameters.threshold=1e-5;
            operators(1).parameters.warmstart=1;
            operators(2).type='sparse_unit-norm';
            operators(2).parameters.threshold=lambdasparseunitnorm;
            operators(3).type='positive';
            optionsadmm.nbitermax=500;
            optionsadmm.seuil=1e-5;
            options.matapprox.leftregparam.optionsadmm=optionsadmm;
            options.matapprox.leftregparam.operatorsadmm=operators;
            toyBenchMSE(comment,learning,options);

    end;
end;

