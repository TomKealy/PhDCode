%  SparsePath -- initialize path to include SparseLab 
%
    global SPARSELABVERSION
    SPARSELABVERSION = 100;
    BVstr=num2str(SPARSELABVERSION);
    
	fprintf('\nWelcome to SparseLab v %g\n\n', SPARSELABVERSION);
%
	global SPARSELABPATH
	global PATHNAMESEPARATOR
	global PREFERIMAGEGRAPHICS
	global MATLABPATHSEPARATOR
	global WLVERBOSE
   
	WLVERBOSE = 'Yes';
	PREFERIMAGEGRAPHICS = 1;
%
	Friend = computer;
	if strcmp(Friend,'MAC2'),
	  PATHNAMESEPARATOR = ':';
	  SPARSELABPATH = ['Macintosh HD:Build:Sparselab', BVstr, PATHNAMESEPARATOR];
	  MATLABPATHSEPARATOR = ';';
	elseif isunix,
	  PATHNAMESEPARATOR = '/';
	  SPARSELABPATH = [pwd, PATHNAMESEPARATOR];
	  MATLABPATHSEPARATOR = ':';
	elseif strcmp(Friend(1:2),'PC');
	  PATHNAMESEPARATOR = '\';	  
	  SPARSELABPATH = [matlabroot,'\toolbox\SparseLab', BVstr, PATHNAMESEPARATOR];  
   	  %SPARSELABPATH = [pwd, PATHNAMESEPARATOR];
	  MATLABPATHSEPARATOR = ';';
	else
		disp('I don''t recognize this computer; ')
		disp('Pathnames not set; solution: edit SparsePath.m\n\n')
	end
%
	global MATLABVERSION
	V = version;
	MATLABVERSION = str2num(V(1:3));

    if MATLABVERSION < 5.3,
        disp('Warning: This version is only supported on Matlab 7.x');
        Sparsep=genpath(SPARSELABPATH,1);
    else
        Sparsep=genpath(SPARSELABPATH);
    end
    
    addpath(Sparsep);   
%
	fprintf('Setting Global Variables:\n');
	fprintf('   global MATLABVERSION = %g\n',            MATLABVERSION)
    fprintf('   global SparseLABVERSION = %g\n',          SPARSELABVERSION)
	fprintf('   global SparseLABPATH = %s\n',                SPARSELABPATH)
	fprintf('   global PATHNAMESEPARATOR = "%s"\n',  PATHNAMESEPARATOR)
	fprintf('   global PREFERIMAGEGRAPHICS = %g\n',PREFERIMAGEGRAPHICS)

    % check if SparseLab is installed successfully
    dr=dir([matlabroot,PATHNAMESEPARATOR,'toolbox']);
    
    l=length(SPARSELABPATH);
    %dr=dir(SPARSELABPATH(1:(l-13)))
    doneset = 0;
    
for k=1:length(dr)
    if (dr(k).name(1)~='.')
        if dr(k).isdir==1
            tmp=min(length(dr(k).name), 12);
            if strcmp(['sparselab', BVstr], lower(dr(k).name(1:tmp)))
                doneset = 1;
	            fprintf('\nSparseLab %g Setup Complete\n\n', SPARSELABVERSION)
                break;
            end
        end
    end
end

if doneset~=1,
    fprintf('\nError occurs and SparseLab %g has not been set up.\n', SPARSELABVERSION)
    fprintf('Solution: edit SparsePath.m\n\n')
    clear all;
    break;
end

	fprintf('Currently available browsers for reproducing figures from the following papers:\n');
%   fprintf('   BrowserAllPaper - Interface to view figures for all papers\n');
    fprintf('   ExtCSDemo - demo for paper "Extensions of Compressed Sensing"\n');
    fprintf('   HDCPNPDDemo - demo for paper "High-Dimensional Centrosymmetric Polytopes with Neighborliness Proportional to Dimension"\n');
    fprintf('   MSNVENODemo - demo for paper "Model Selection When the Number of Variables Exceeds the Number of Observations"\n');
    fprintf('   NPSSULEDemo - demo for paper "Neighborly Polytopes and Sparse Solutions of Underdetermined Linear Equations"\n');
    fprintf('   NRPSHDDemo - demo for paper "Neighborliness of Randomly-Projected Simplices in High Dimensions"\n');
    fprintf('   SNSULELPDemo - demo for paper "Sparse Nonnegative Solutions of Underdetermined Linear Equations by Linear Programming"\n');
    fprintf('   StOMPDemo - demo for paper "Sparse Solution of Underdetermined Linear Equations by Stagewise Orthogonal Matching Pursuit"\n'); 
    
    fprintf('Currently available examples:\n');
    fprintf('   Nonnegative Factorization\n');
    fprintf('   Signal Reconstruction\n');
    fprintf('   Regression Example\n');
    fprintf('   Time-Frequency Separation\n');
    
    fprintf('\nFor more information, please visit: \n');
    fprintf('   http://sparselab.stanford.edu\n\n');
    
%if strcmp(WAVELABINSTALLED, 'NO')
    fprintf('Please ignore the following message if WaveLab has been installed.\n\n');
    fprintf('There are SparseLab functions which call WaveLab functions.\n');
    fprintf('We recommend that the users download WaveLab from the website\n');
    disp('   http://www-stat.stanford.edu/~wavelab');
    fprintf('and install the package in the directory\n');
    disp(['   ',strcat(matlabroot,PATHNAMESEPARATOR,'toolbox')]);
    fprintf('\n\n');
    %end
    
%
% the next statement leaves items in global workspace
%  but hides them from local workspace
%
%	clear SPARSELABPATH MATLABVERSION PATHNAMESEPARATOR
    clear dr V k BVstr tmp Sparsep doneset l
	clear Friend PREFERIMAGEGRAPHICS MATLABPATHSEPARATOR

%
% Copyright (c) 2006. Victoria Stodden and David Donoho
% 

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
