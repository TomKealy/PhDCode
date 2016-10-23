function [ w, J ] = activesetalgorithm( X, XtY_over_n, G, D, IDX, lambda, params)
                                                 
%ACTIVE SET ALGORITHM FOR GENERAL GROUPS ON 2-D GRIDS
%
% [ w, J ] = activesetalgorithm( X, XtY_over_n, G, D, IDX, lambda, params)
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
% params      : structure that contains additional paramaters (see the fields below)
%
% OUTPUTS:
% 
% w           : loading vector
% J           : support


[n,p] = size(X);

nd = size(IDX,1);% number of group directions, for instance nd = 4 for rectangles

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
if isfield( params, 'constrained_mode' ),% Determines wether the constrained formulation is used
                                          % By default (as in the paper),
                                          % the squared penalization is used 
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
            sufficient_condition_test = sqrt( lambda*2*epsilon );
        
        else 
            
            necessary_condition_test  = 0;
            sufficient_condition_test = epsilon / lambda;
            
        end
        
        
    else
        
        grad_dot_w = grad(J)'*w(J);
        
        if opt_params.squared_penalization_mode,
            
            necessary_condition_test  = sqrt( -lambda*grad_dot_w );
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
            
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            
            beyond_fringe_group_idx = smallest_GJ_idx(  smallest_GJ_idx+1 < IDX(:,2)  ) + 2;
            
            free_directions         = find( smallest_GJ_idx < IDX(:,2) );
            
            fringe_group_idx        = smallest_GJ_idx( free_directions ) + 1;

            fringe_variable_idx     = find( ~( any( G(:, beyond_fringe_group_idx ), 2 ) | J_bool ) );
            
            reduced_G               = G( fringe_variable_idx, fringe_group_idx) ;
            

            % K is length(fringe_variable_idx) x n_patterns logical matrix
            K                       = mexGetFringeNonzeroPatterns( reduced_G );
            
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            
            norm_grad_over_sum_dG = zeros( size(K,2), 1);
            G_K                   = false( length(free_directions), size(K,2) );
            
            for j=1:size(K,2),
                
               tmp1     = fringe_variable_idx( K( :, j) );
               
               G_K(:,j) = any( reduced_G(  K( :, j), : ), 1 )';
               
               tmp2     = fringe_group_idx( G_K( :, j) );

               % We calculate the max of the grad-norm-over-sum-norminf dG 
               sum_norm_inf_dG          = sum(  max(  D( tmp1, tmp2 ), [], 1 )  ); % max over rows
               
               norm_grad_over_sum_dG(j) = norm( grad( tmp1 ), 2 ) / sum_norm_inf_dG;

            end;
            
            [max_, argmax_in_fringe_group_idx ] =  max( norm_grad_over_sum_dG );

            J_to_add = fringe_variable_idx( K(:,argmax_in_fringe_group_idx) );
             
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
                
                directions_to_update = free_directions( G_K( :,argmax_in_fringe_group_idx) );
                
                smallest_GJ_idx( directions_to_update ) = smallest_GJ_idx( directions_to_update ) + 1;
                
            end
            
            
            J                 = [J; J_to_add]; %#ok<AGROW>
            J_bool(J_to_add)  = true;
            
        end
        
        
        
    else

        %==================================================================
        % Sufficient condition 
        %==================================================================

        beyond_fringe_group_idx = smallest_GJ_idx(  smallest_GJ_idx+1 < IDX(:,2)  ) + 2;

        free_directions         = find( smallest_GJ_idx < IDX(:,2) );

        fringe_group_idx        = smallest_GJ_idx( free_directions ) + 1;

        fringe_variable_idx     = find( ~( any( G(:, beyond_fringe_group_idx ), 2 ) | J_bool ) );

        reduced_G               = G( fringe_variable_idx, fringe_group_idx) ;
        
        % K is length(fringe_variable_idx) x n_patterns logical matrix
        
        K                       = mexGetFringeNonzeroPatterns( reduced_G );
        
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        
        sum_dG = sum( D( :, mexMultiColon( smallest_GJ_idx+1, IDX(:,2) ) ), 2);
        
        bound  = zeros( length(fringe_group_idx), 1);
        
        for j=1:length(fringe_group_idx),

            g        = G( :, fringe_group_idx(j) );
            
            bound(j) = norm( grad( g ) ./ sum_dG( g ),  2);
               
        end   
        
        %------------------------------------------------------------------
        
        [max_, argmax_in_fringe_group_idx ] =  max( bound );
        
        G_star = G( fringe_variable_idx, fringe_group_idx(argmax_in_fringe_group_idx) );
        
        % What are the patterns in K that intersect G_star ?
        fringe_pattern_inter_G_star = any( mexAndLogic( G_star, K ), 1 );
        
        
        if ~any(fringe_pattern_inter_G_star),
            
            % We look at all the groups that intersect G_star and then
            % apply the same scheme as before
            reduced_G_inter_G_star = any( mexAndLogic( G_star, reduced_G ), 1 );
            
            for g=find(reduced_G_inter_G_star),
                fringe_pattern_inter_G_star = fringe_pattern_inter_G_star | any( mexAndLogic( reduced_G(:,g), K ), 1 );
            end
            
        end
        
        union_of_patterns_in_K = any( K( :, fringe_pattern_inter_G_star ), 2 );
        J_to_add               = fringe_variable_idx( union_of_patterns_in_K );
        %------------------------------------------------------------------
        %------------------------------------------------------------------       
        
        if (  max_ <= sufficient_condition_test ),        
            SufficientCondition = true %#ok<NOPRT>
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
                
                G_K                                     = any( reduced_G( union_of_patterns_in_K, : ), 1 );
                
                directions_to_update                    = free_directions( G_K );
                
                smallest_GJ_idx( directions_to_update ) = smallest_GJ_idx( directions_to_update ) + 1;
                
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
    
    
    w(J) = slasso_SDPT3( XtX_J_over_n(J,:), XtY_over_n(J), G( J, active_G_idx ), D( J, active_G_idx ), lambda, opt_params);
    
    %======================================================================
    
    t = t+1;
   
end

end

%==========================================================================
%==========================================================================
