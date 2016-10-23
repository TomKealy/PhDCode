function [ w, J ] = activesetalgorithm_rectangle( X, XtY_over_n, G, D, IDX, lambda, Height, Width, params)
                                                 
%ACTIVE SET ALGORITHM FOR RECTANGULAR GROUPS ON 2-D GRIDS
%
% [ w, J ] = activesetalgorithm_rectangle( X, XtY_over_n, G, D, IDX, lambda, Height, Width, params)
%
% INPUTS:
%
% X           : design matrix (size n x p, i.e., n row data of dimension p)
% XtY_over_n  : X'*Y/n, correlations between design matrix and output vector Y (size p x 1)
% G           : group (logical) matrix (size p x ng, with ng the number of groups)
%               G(i,j) = true if variable i is in group j, false otherwise
% D           : weight matrix (same size as G)
%               D(i,j) = weight of the variable i in the group j if i is in j, 0 otherwise
% IDX         : index matrix that contains the sorting of G, per direction and size
%               (size nd x 2, with nd the number of directions in G; for instance 4 for the rectangular groups)
%               More details in the description of the function get_groups
% lambda      : regularization parameter
% Height      : height of the grid
% Width       : width of the grid
% params      : structure that contains additional paramaters (see the fields below)
%
% OUTPUTS:
% 
% w           : loading vector
% J           : support
%
%   Copyright (c) 2009 Rodolphe Jenatton. All rights reserved.




[n,p] = size(X);

nd = size(IDX,1);% number of group directions, here nd = 4

%==========================================================================
% Process params
%--------------------------------------------------------------------------
if isfield( params, 'maxCardJ' ),% Max size of the support
    maxCardJ = params.maxCardJ;
else
    maxCardJ = floor(p/10);
end
%--------------------------------------------------------------------------
if isfield( params, 'epsilon' ), % Precision of the optimization
    epsilon = params.epsilon;
else
    epsilon = 1e-3;
end
%--------------------------------------------------------------------------
if isfield( params, 'MaxIterations' ), % Max number of iterations
    MaxIterations = params.MaxIterations;
else
    MaxIterations = 15;
end
%--------------------------------------------------------------------------
if isfield( params, 'constrained_mode' ), % Determines wether the constrained formulation is used
                                          % By default (as in the paper), the squared penalization is used 
    opt_params.constrained_mode          = params.constrained_mode;
    opt_params.squared_penalization_mode = ~opt_params.constrained_mode;
else
    opt_params.squared_penalization_mode = true;
    opt_params.constrained_mode          = false;
end
%==========================================================================
% Initialization

J               = [];
J_bool          = false(p,1);
smallest_GJ_idx = zeros(nd,1);

w = zeros(p,1);

grad_over_sum_dG_norm = zeros(nd,1);

t = 1;

NecessaryCondition  = false;
SufficientCondition = false;

%==========================================================================

while ~( NecessaryCondition && SufficientCondition ) && ... 
       ( t <= MaxIterations )                        && ...
       ( length(J) < maxCardJ ),
    
    %======================================================================
    % Warning, it is minus the gradient of the loss function  
    
    if isempty(J),
        grad = -XtY_over_n;
    else
        grad = -XtY_over_n + XtX_J_over_n*w(J); 
    end;
    
    %======================================================================
    % Calculation of the quantity used in the tests for the NecessaryCondition/SufficientCondition
    
    if isempty(J),
        
        if opt_params.squared_penalization_mode,
            
            necessary_condition_test  = 0;
            sufficient_condition_test = sqrt( lambda*( 2*epsilon ) );
        
        else 
            
            necessary_condition_test  = 0;
            sufficient_condition_test = epsilon / lambda;
            
        end
        
        
    else
        
        grad_dot_w = grad(J)'*w(J);
        
        if opt_params.squared_penalization_mode,
            
            necessary_condition_test  = sqrt( -lambda * grad_dot_w );
            sufficient_condition_test = sqrt( lambda*( 2*epsilon -grad_dot_w ) );
        
        else 
            
            necessary_condition_test  = -grad_dot_w / lambda;
            sufficient_condition_test = (epsilon-grad_dot_w) / lambda;
            
        end
        
        
    end
    %======================================================================
    % Necessary condition 
    %======================================================================
    
    if ~NecessaryCondition, 
    
    
        if isempty(J),
            
            % We need to find the leaves of the DAG of nonzero patterns
            % ---- WE ASSUME THE LEAVES ARE THE SINGLETONS ----
          
            [max_, J_to_add ] = max( abs(grad) ./ sum( D, 2 ) );          
            
        else
            
            %   ---- WARNING: convention ----
            %
            %   G{1} : horizontal groups from right to left
            %   G{2} : horizontal groups from left to right
            %   G{3} : vertical   groups from bottom to top
            %   G{4} : vertical   groups from top to bottom
            
            fringe_variable_idx = get_fringe_variables( J, Height, Width ); %fringe_variable_idx is p x 4 logical matrix
             
            for d=1:nd,
                
                K = fringe_variable_idx(:,d);
                
                if any( K ),
                    % We calculate the max of the grad-norm-over-sum-norminf-dG 
                    
                    grad_over_sum_dG_norm(d) = norm( grad(K), 2) / max( D( K, smallest_GJ_idx(d)+1 ), [], 1 );% Max over rows
                else
                    % We put -1 by default not to achieve the max
                    grad_over_sum_dG_norm(d) = -1;  
                end
                                   
            end;

            [max_, argmax_in_fringe_group_idx ] =  max( grad_over_sum_dG_norm );

            J_to_add =  find( fringe_variable_idx(:,argmax_in_fringe_group_idx) );
             
        end
        %------------------------------------------------------------------
        %------------------------------------------------------------------
   
        if (  max_ <= necessary_condition_test ),
            NecessaryCondition = true %#ok<NOPRT>
        else
            
            if isempty(J),
                
                for d=1:nd,
                    
                    shift_idx = find( any( G( J_to_add, IDX(d,1):IDX(d,2) ), 1 ), 1, 'last');
                    
                    if isempty(shift_idx),
                        smallest_GJ_idx(d) = IDX(d,1); 
                    else
                        smallest_GJ_idx(d) = IDX(d,1) +shift_idx- 1;
                    end
                    
                end
                
            else
                
                smallest_GJ_idx(argmax_in_fringe_group_idx) = ...
                    smallest_GJ_idx(argmax_in_fringe_group_idx) + 1;
                
            end
            
            
            J                 = [J; J_to_add];%#ok<AGROW,NOPRT>
            J_bool(J_to_add)  = true;
            
        end
        
        
        
    else

        %==================================================================
        % Sufficient condition 
        %==================================================================
   
        sum_dG = sum( D( :, mexMultiColon( smallest_GJ_idx+1, IDX(:,2) ) ) , 2);
        
        for d=1:nd,
            
            if smallest_GJ_idx(d) < IDX(d,2),
            
                g  = G( :, smallest_GJ_idx(d)+1 );% variables in G, for G in Fringe.
                  
                grad_over_sum_dG_norm(d) = norm( grad( g ) ./ sum_dG( g ),  2);
                
            else
                
                grad_over_sum_dG_norm(d) = -1;
                
            end
            
        end    

        [max_, argmax_in_fringe_group_idx ] = max( grad_over_sum_dG_norm );
        
 
        fringe_variable_idx = get_fringe_variables( J, Height, Width ); %fringe_variable_idx is p x 4
        
 
        J_to_add =  find( fringe_variable_idx( :, argmax_in_fringe_group_idx ) );
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        
        if (  max_ <= sufficient_condition_test ),        
            SufficientCondition = true %#ok<NOPRT>
        else
            
            if isempty(J),
                
                shift_idx = find( any( G( J_to_add, IDX(d,1):IDX(d,2) ), 1 ), 1, 'last');

                if isempty(shift_idx),
                    smallest_GJ_idx(d) = IDX(d,1); 
                else
                    smallest_GJ_idx(d) = IDX(d,1) +shift_idx- 1;
                end
                
            else
                
                smallest_GJ_idx(argmax_in_fringe_group_idx) = ...
                    smallest_GJ_idx(argmax_in_fringe_group_idx) + 1;
                
            end
            
            J                 = [J; J_to_add]; %#ok<AGROW,NOPRT>
            J_bool(J_to_add)  = true;
            
        end
        %==================================================================
        
    end
   
    %======================================================================
    % Calculation of the covariance matrix reduced to J 
    
    XtX_J_over_n = X'*X(:,J)/n;
    
    %======================================================================
    % Solve reduced problem on J
     
    active_G_idx = mexMultiColon( IDX(:,1), smallest_GJ_idx );
       
    w(J) = slasso_SDPT3(XtX_J_over_n(J,:), XtY_over_n(J), G(J,active_G_idx), D(J,active_G_idx), lambda, opt_params);
    
    %======================================================================
    
    t = t+1;
   
end

end

%==========================================================================
%==========================================================================

function K = get_fringe_variables( J, H, L )

% J is the active set
% H is the height of the grid
% L is the width of the grid
    
    p = H * L;
    
    K = false(p,4);
    
    Card_J = length(J);

    min_J = min(J);
    max_J = max(J);
    
    % l : width of J
    % h : height of J
    
    l = ( H + 1 + max_J-min_J + sqrt( (H + 1 + max_J-min_J)^2 - 4*H*Card_J ) ) / (2*H);
    h = Card_J / l;
   
    
    if  max_J <= (L-1)*H,% there is a column on the right
       
        K( ( min_J+l*H ):( min_J+l*H + h-1 ), 1 ) = true;
  
    end 
    

    if  min_J > H, % there is a column on the left
       
        K( (min_J-H):(min_J + h-1 -H), 2 ) = true;
         
    end

    
    if  mod(min_J + h-1, H) > 0, % there is a row below
       
        K( (min_J + h):H:( min_J + h + (l-1)*H ), 3 ) = true;
               
    end
       

    if  mod(min_J,H) > 1, % there is a row above
       
        K( (min_J - 1):H:(min_J - 1 + (l-1)*H ), 4 ) = true;
               
    end
 
 
end
