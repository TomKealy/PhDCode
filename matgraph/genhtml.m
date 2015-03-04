function genhtml
% generate matgraph html documentation using m2html
tic
addpath('/Users/ers/Common/Programming/matlab/m2html/')
save_dir = pwd;
cd ~/Common/Programming/

file_list = cell(1,4);
file_list{1} = 'matgraph';
file_list{2} = 'matgraph/@graph';
file_list{3} = 'matgraph/@partition';
file_list{4} = 'matgraph/@permutation';

m2html('mfiles',file_list, ...
    'htmldir','matgraph/html', ...
    'recursive','off' ...
    )

cd matgraph
rmpath('/Users/ers/Common/Programming/matlab/m2html/')
et = toc;
mins = floor(et/60);
secs = round(et - 60*mins);

disp(['Elapsed time: ', int2str(mins), '''', int2str(secs),'"'])
cd(save_dir)