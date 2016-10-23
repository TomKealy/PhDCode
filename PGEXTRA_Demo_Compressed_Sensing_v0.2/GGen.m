function [kappa_G,P,per_t,fig]=GGen(L,per,name,show_graph,quiet,arg_opt1)
% -------------------------------------------------------------------------
%	GGen randomly generates certain class of graphs with specified connectivity
%	ratios.
%
%	Required parameters:
%       L: number of nodes
%       per: expected connectivity ratio
%       name: the type of graph, can be random, Line, Star, Complete, Cycle, Regular, Bipartite, 
%				Biclique, k-NNG, RGG, Hierarchical (Hierarchical has not been implemented yet),
%				2D grid, and 3D grid
%       show_graph: plot the graph or not, true/false
%       quiet: silent mode or not, true/false
%       arg_opt1: used for specific classes of random graphs
%           (1) Regular graph: arg_opt1 is the degree of the nodes
%           (2) Bipartite graph: the number of nodes in one part
%           (3) k-NNG graph: k, the number neighbours
%           (4) RGG: the Radius of connectivity
%           (5) 2D Grid: a vector containing the length and the width of the 2D Grid
%           (6) 3D Grid: a vector containing the length, the width, and the height of the 3D Grid
%   
%   Return values:
%       kappa_G: the condition number of the graph
%       P: Adjacency matrix of the generated graph
%       per_t: the actual connectivity ratio
%       fig: the figure handle of the plotted graph, exists only when show_graph=true
%
%	Examples:
%       [kappa_G,P,per_t,fig]=GGen(10,0.5,'random',true,true);
%
%   References: 
%	[1] W. Shi, Q. Ling, K. Yuan, G. Wu, and W. Yin, ¡°On the Linear
%	Convergence of the ADMM in Decentralized Consensus Optimization,¡±
%	Signal Processing, IEEE Transactions on, vol. 62, no. 7, pp. 1750-1761, 2014.

%   Copyright, Wei Shi, Institute of Industry Automation, University of
%   Science and Technology of China
%   Please use the code under the GNU GPLv3
% ------------------------------------------------------------------------- 
L_A=0;
L_B=0;
P(1:L,1:L)=0;
if strcmp(name,'Regular')
    if ~quiet
        warning('Generating random regular graph needs matgraph! By now, we can only generate random regular graph with degree at most 6 efficiently!');
    end;
    d_for_regular=arg_opt1;
elseif strcmp(name,'Bipartite')||strcmp(name,'Complete Bipartite')||strcmp(name,'Biclique')
    L_A=arg_opt1;
elseif strcmp(name,'Random Proximity')||strcmp(name,'k-NNG')
    Neighbours=arg_opt1;
elseif strcmp(name,'Random Geometric')||strcmp(name,'RGG')
    Radius=arg_opt1;
elseif strcmp(name,'Hierarchical')
    L_Star=arg_opt1; %#ok<NASGU>
elseif strcmp(name,'2D Grid')
    m=arg_opt1(1);
    n=arg_opt1(2);
elseif strcmp(name,'3D Grid')
    m=arg_opt1(1);
    n=arg_opt1(2);
    p=arg_opt1(3);
end;

if (~exist('name','var'))||(strcmp(name,'random'))   
    if per<2/L   
        selected_set_length=0;
        per_t=2/L;
        if ~quiet
            fprintf('Ratio has been reset as %f .\n',per_t);
        end
    elseif per>1 
        selected_set_length=L*(L-1)/2-L+1;
        per_t=1;
        if ~quiet
            fprintf('Ratio has been reset as %f .\n',per_t);
        end;
    else
        selected_set_length=round(per*L*(L-1)/2-L+1);
        per_t=2*(selected_set_length+L-1)/(L*(L-1));
        if ~quiet
            fprintf('Ratio is %f .\n',per_t);
        end;
    end;
    vertices=1:L;
    for i=L:-1:2
        temp_edge=vertices(randperm(i,2));
        P(temp_edge(1),temp_edge(2))=1;
        P(temp_edge(2),temp_edge(1))=1;
%         vertices=setdiff(vertices,temp_edge(1));
        for temp_i=1:i
            if temp_edge(1)==vertices(temp_i)
                vertices=[vertices(1:temp_i-1),vertices(temp_i+1:end)];
                break;
            end;
        end;
    end;
    temp_A=P+tril(~P);
    candidate_set=find(~temp_A);
    selected_set=candidate_set(randperm(length(candidate_set),selected_set_length));
    temp_A(selected_set)=1;
    P=triu(temp_A,1);
    P=P+P';
elseif strcmp(name,'Line')
    if ~quiet
        fprintf('A "Line" generated.\n');
    end;
    per_t=2/L;
    for tempi=1:(L-1)
        P(tempi,tempi+1)=1;
    end;
    for tempi=2:L
        P(tempi,tempi-1)=1;
    end;
elseif strcmp(name,'Star')
    if ~quiet
        fprintf('A "Star" generated.\n');    
    end;
    per_t=2/L;
    tempi=1;
    for tempj=2:L
        P(tempi,tempj)=1;
    end;
    tempj=1;
    for tempi=2:L
        P(tempi,tempj)=1;
    end;
elseif strcmp(name,'Complete')
    if ~quiet
        fprintf('A complete graph generated.\n');  
    end;
    per_t=1;
    P=ones(L,L)-diag(ones(L,1));
elseif strcmp(name,'Cycle')
    if ~quiet
        fprintf('A "Cycle" generated.\n');
    end;
    per_t=2/(L-1);
    tempcount=0;
    for tempi=1:(L-1)   
        P(tempi,tempi+1)=1;
    end;
    P(L,1)=1;
    for tempi=2:L
        tempcount=tempcount+1;
        P(tempi,tempi-1)=1;
    end; 
    P(1,L)=1;
elseif strcmp(name,'Regular')
    addpath('matgraph');
    % graph_destroy;
    graph_init(1000);
    if ~quiet
        fprintf('A regular graph generated.\n');
    end;
    if d_for_regular>min(L-1,6);
        d_for_regular=min(L-1,6);
        if ~quiet
            fprintf(['Degree has been reset as ',num2str(d_for_regular),' .\n']);
        end;
    elseif d_for_regular<2
        d_for_regular=2;
        if ~quiet
            fprintf('Degree has been reset as 2 .\n');
        end;
    elseif mod(L*d_for_regular,2)==1
        d_for_regular=d_for_regular+sign(random('unif',-1,1));
        if ~quiet
            fprintf(['Degree has been reset as ',num2str(d_for_regular),' .\n']);
        end;
    else
        if ~quiet
            fprintf(['Degree is ',num2str(d_for_regular),' .\n']);
        end;
    end;
    per_t=d_for_regular/(L-1);
    g=graph;
    connected_flag=false;
    while ~connected_flag
        random_regular(g,L,d_for_regular);
        connected_flag=isconnected(g);
    end;
    P=matrix(g);
elseif strcmp(name,'Bipartite')
    if ~quiet
        fprintf('A bipartite graph generated.\n');
    end;
    if L<3
        if ~quiet
            warning('Nodes number is smaller than 3!');
        end;
    else
        L_B=L-L_A;
        if L_B<L_A
            L_A=L_B;
            L_B=L-L_A;
        end;
    end;
    if per>2*L_A*L_B/(L*(L-1))
        per_t=2*L_A*L_B/(L*(L-1));
        edges_num=L_A*L_B;
        if ~quiet
            fprintf('Ratio has been reset as %f .\n',per_t);
        end;
    elseif per<2/L
        per_t=2/L;
        edges_num=L-1;
        if ~quiet
            fprintf('Ratio has been reset as %f .\n',per_t);
        end;
    else
        edges_num=round(per*L*(L-1)/2);
        per_t=2*edges_num/(L*(L-1));
        if ~quiet
            fprintf('Ratio is %f .\n',per_t);
        end;
    end;
    P(1:L,1:L)=0;
    Set_A=1;
    len_A=1;
    candidate_A=2:L_A;
    len_can_A=L_A-1;
    Set_B=L_A+1;
    len_B=1;
    candidate_B=L_A+2:L;
    len_can_B=L_B-1;
    P(1,L_A+1)=1;
    P(L_A+1,1)=1;
    edges_num=edges_num-1;
    temp_edge(1:2)=0;
    for i=1:L-2
        if len_can_A/(len_can_A+len_can_B)>rand
            temp_edge(1)=candidate_A(randperm(len_can_A,1));
            temp_edge(2)=Set_B(randperm(len_B,1));
            P(temp_edge(1),temp_edge(2))=1;
            P(temp_edge(2),temp_edge(1))=1;
            Set_A=[Set_A,temp_edge(1)]; %#ok<AGROW>
            len_A=len_A+1;
            candidate_A=setdiff(candidate_A,temp_edge(1));
            len_can_A=len_can_A-1;
            edges_num=edges_num-1;
        else
            temp_edge(1)=candidate_B(randperm(len_can_B,1));
            temp_edge(2)=Set_A(randperm(len_A,1));
            P(temp_edge(1),temp_edge(2))=1;
            P(temp_edge(2),temp_edge(1))=1;
            Set_B=[Set_B,temp_edge(1)]; %#ok<AGROW>
            len_B=len_B+1;
            candidate_B=setdiff(candidate_B,temp_edge(1));
            len_can_B=len_can_B-1;
            edges_num=edges_num-1;
        end;
    end;
    if edges_num>0
        temp_P=P(1:L_A,L_A+1:L);
        candidate_set=find(temp_P==0);
        temp_P(candidate_set(randperm(length(candidate_set),edges_num)))=1;
        P=[zeros(L_A),temp_P;temp_P',zeros(L_B)];
    end;
elseif strcmp(name,'Complete Bipartite')||strcmp(name,'Biclique')
    if ~quiet
        fprintf('A complete bipartite graph generated.\n');
    end;
    L_B=L-L_A;
    per_t=2*L_A*L_B/(L*(L-1));
    Set_A=1:L_A;
    Set_B=(L_A+1):L;
    P(Set_A,Set_B)=1;
    P(Set_B,Set_A)=1;
elseif strcmp(name,'Random Proximity')||strcmp(name,'k-NNG')   
    if Neighbours>L-1
        Neighbours=L-1;
        if ~quiet
            fprintf(['The K-NNG has been reset as ',num2str(Neighbours), '-NNG.\n']);
        end;
    elseif Neighbours<1
        Neighbours=1;
        if ~quiet
            fprintf(['The K-NNG has been reset as ',num2str(Neighbours), '-NNG.\n']);
        end;
    else
        if ~quiet
            fprintf(['The K-NNG is ',num2str(Neighbours), '-NNG.\n']);
        end;
    end;
    V=random('unif',0,1,L,3);
    EDM=squareform(pdist(V,'euclidean'));
    [~,iEDM]=sort(EDM,'ascend');
    jEDM=kron((iEDM(1,:))',ones(Neighbours,1));
    iEDM=reshape(iEDM(2:Neighbours+1,:),[],1);   
    P((iEDM-1)*L+jEDM)=1;
    P((jEDM-1)*L+iEDM)=1;
    per_t=sum(P(:))/(L*(L-1));
    fprintf('per is %f .\n',per_t);
    g=graph;
    fast_set_matrix(g,P);
    if isconnected(g)==0
        if ~quiet
            warning('Graph not connected!');
            per_t=-1;
        end;
    end;  
elseif strcmp(name,'Random Geometric')||strcmp(name,'RGG')
    temp_L=L;
    V=[];
    while temp_L>0
        temp_V=random('unif',-1,1,ceil(temp_L*6/pi),3);
        temp_V=temp_V(sum(temp_V.^2,2)<=1,:);
        temp_L=L-size(V,1)-size(temp_V,1);
        if temp_L>=0
            V=[V;temp_V]; %#ok<AGROW>
        else
            V=[V;temp_V(1:L-size(V,1),:)]; %#ok<AGROW>
        end;
    end;
%     figure,plot3(V(:,1),V(:,2),V(:,3),'o');%for debugging;
    EDM=squareform(pdist(V,'euclidean'))+diag(ones(L,1)*inf);
    P(EDM<Radius)=1;   
    per_t=sum(P(:))/(L*(L-1));
    fprintf('per is %f .\n',per_t);
    g=graph;
    fast_set_matrix(g,P);
    if isconnected(g)==0
        if ~quiet
            warning('Graph not connected!');
            per_t=-1;
        end;
    end;
elseif strcmp(name,'2D Grid')
    L=m*n;
    P=zeros(L,L);
    G_2D=reshape(randperm(L),m,n);    
    V1=G_2D(1:m-1,:);
    V2=G_2D(2:m,:);
    len_V=length(V1(:));
    for k=1:len_V
        P(V1(k),V2(k))=1;
        P(V2(k),V1(k))=1;
    end;
    V1=G_2D(:,1:n-1);
    V2=G_2D(:,2:n);
    len_V=length(V1(:));
    for k=1:len_V
        P(V1(k),V2(k))=1;
        P(V2(k),V1(k))=1;
    end;
    per_t=2*((m-1)*n+m*(n-1))/(L*(L-1));
    fprintf('per is %f .\n',per_t);
elseif strcmp(name,'3D Grid')
    L=m*n*p;
    P=zeros(L,L);
    G_3D=reshape(randperm(L),m,n,p);    
    V1=G_3D(1:m-1,:,:);
    V2=G_3D(2:m,:,:);
    len_V=length(V1(:));
    for k=1:len_V
        P(V1(k),V2(k))=1;
        P(V2(k),V1(k))=1;
    end;
    V1=G_3D(:,1:n-1,:);
    V2=G_3D(:,2:n,:);
    len_V=length(V1(:));
    for k=1:len_V
        P(V1(k),V2(k))=1;
        P(V2(k),V1(k))=1;
    end;
    V1=G_3D(:,:,1:p-1);
    V2=G_3D(:,:,2:p);
    len_V=length(V1(:));
    for k=1:len_V
        P(V1(k),V2(k))=1;
        P(V2(k),V1(k))=1;
    end;
    per_t=2*((m-1)*n*p+m*(n-1)*p+m*n*(p-1))/(L*(L-1));
    fprintf('per is %f .\n',per_t);
end;

D=diag(sum(P));
L_neg=D-P;
L_pos=abs(L_neg);
svd_L_pos=svd(L_pos);
svd_L_neg=svd(L_neg);
kappa_G=(svd_L_pos(1)/svd_L_neg(end-1))^(1/2);
if (exist('show_graph','var'))&&(show_graph) 
    ids={};
    for i=1:L
        ids=union(ids,num2str(i));
    end;
    if strcmp(name,'Bipartite')||strcmp(name,'Complete Bipartite')||strcmp(name,'Biclique')
        bg=biograph(tril(P),ids,'ShowArrows','off','LayoutType','hierarchical','EdgeType','curved');
    else
        bg=biograph(tril(P),ids,'ShowArrows','off','LayoutType','equilibrium','EdgeType','curved');   
    end;
    Nodes=get(bg,'Nodes');
    set(Nodes,'Shape','circle','LineWidth',2,'Size',[14,14],'FontSize',12);  
    Edges=get(bg,'Edges');
    set(Edges,'LineColor',[0 0 0],'LineWidth',2)
    bg_gui=biograph.bggui(bg);
    fig=get(bg_gui.biograph.hgAxes,'parent');
else
    fig=0;
end;
end