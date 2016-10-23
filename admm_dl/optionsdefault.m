
addpath('../ksvd/');
addpath('../ksvd/');
addpath('../ksvd/ompbox/');
addpath('../ksvd/tools/');
addpath('..')

%--------------------------------------------------------------------------
%           Problem parameters
%--------------------------------------------------------------------------
learning.dim=50;        % problem dimension
learning.TailleDic=100;   % Dictionary size
learning.T=3;            % nb de fonctions de bases
learning.K=200;           % nb de signaux
learning.nbitermax=10;   
learning.rsnr=5;
learning.verbose=1;
learning.type='smooth'
learning.nbtest=5000;
learning.nbval=100;
learning.lambdasparsecodevec=[0.4]; % regularization parameters
learning.lambdadico=0;   
learning.lambdasparsecodevectest=[0.1];
comment='';

%--------------------------------------------------------------------------
%       Dictionary Learning Options
%--------------------------------------------------------------------------

options.nbitermax=100;                   % nb of iterations 
options.threshold=1e-4;                 % threshold on relative variation 
options.warmstart=0;        
options.eta=1;                          % multiplicative decrease of stopping 
                                        % criterion
options.ALnbitermax=1000;
options.ALmu=10;
options.ALmudico=1;
%--------------------------------------------------------------------------                                         
%  Majorization method options (Iterative Thresholding algorithm)
%--------------------------------------------------------------------------
options.it.threshold=1e-6;
options.it.nbitermax=200000;

%--------------------------------------------------------------------------
% Fista / Augmented Lagrangian options  
%--------------------------------------------------------------------------
% matrix approximation options
options.matapprox.nbitermax=1000;
options.matapprox.threshold=1e-4;  % avec tv faut prendre + petit
options.matapprox.thresholdmin=1e-4;

% regularizer options 
options.matapprox.rightregularizer='\ell_1-\ell_1';
options.matapprox.rightregparam=0;

%-------------------------------------------------------------------------
%           dictionary options
%--------------------------------------------------------------------------
% options.matapprox.leftregularizer='proj-unitnorm2';
% options.matapprox.leftregparam=0;


options.matapprox.leftregularizer='admmvect';
% operators(1).type='tv';
% operators(1).parameters.bound=0.1;   % a diminuer si moins de bruit
% operators(1).parameters.nbitertv=300;
% operators(1).parameters.threshold=1e-5;
% operators(1).parameters.warmstart=1;
%operators(2).type='unit-norm';
        operators(1).type='sparse_unit-norm';
         operators(1).parameters.threshold=0.05;
% operators(2).type='positive';
optionsadmm.nbitermax=300;
optionsadmm.seuil=1e-5;
options.matapprox.leftregparam.optionsadmm=optionsadmm;
options.matapprox.leftregparam.operatorsadmm=operators;

