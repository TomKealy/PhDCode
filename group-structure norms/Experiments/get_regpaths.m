function [ prediction_error, estimation_error, hull_estimation_error, cpu_time, output_params ] = get_regpaths( height, width, hull_pattern, Lambda, n_tr, n_runs, additional_arguments )

%==========================================================================
%GET_REGPATHS

% additional_arguments:
%
%   additional_arguments.G
%   additional_arguments.D
%   additional_arguments.IDX
%   additional_arguments.model (either 'lasso', 'slasso', 'islasso', 'activeset_rectangle', 'activeset' or 'slasso_sq_SDPT3' )

%==========================================================================
% Create output_params

output_params.SNR           = 3;
output_params.state         = 0;
output_params.n_tr          = n_tr;
output_params.n_test        = 500;
output_params.n_runs        = n_runs;
output_params.threshold     = 1e-7;
output_params.Kfolds        = 5;
output_params.Lambda        = Lambda;
output_params.hull_pattern  = hull_pattern;
output_params.height        = height;
output_params.width         = width;
output_params.maxCardJ      = 5*length( hull_pattern );
output_params.MaxIterations = 30;

%==========================================================================

randn('state',  output_params.state);
rand('twister', output_params.state);


prediction_error      = zeros(output_params.n_runs,1);
estimation_error      = zeros(output_params.n_runs,1);
hull_estimation_error = zeros(output_params.n_runs,1);
cpu_time              = zeros(output_params.n_runs*length(Lambda)*output_params.Kfolds,1);

%==========================================================================

p =  height * width;

for t=1:output_params.n_runs,
    
    fprintf('\n Run %d out of %d \n', t, output_params.n_runs);
    
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Generate the data: w, X, Y training/test
    
    w        = zeros(p,1);
    
    JHull    = get_random_hull_position( hull_pattern, height, width );
    w(JHull) = randn( length(hull_pattern), 1 );
    
    X_tr  = randn(n_tr,p);
    X_w   = X_tr*w;
    noise = randn(n_tr,1);
    sigma = norm(X_w)/(output_params.SNR*norm(noise));
    Y_tr  = X_w + sigma*noise;
    
    XtX_tr_over_n = X_tr' * X_tr / n_tr;
    XtY_tr_over_n = X_tr' * Y_tr / n_tr;
    
    X_test           = randn(output_params.n_test,p);
    XtX_test_over_n  = X_test' * X_test / output_params.n_test;
    
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Compute lambda_star by CV
    
    [lambda_star, cpu_time_tmp] = ...
        get_lambda_by_CV( w, X_tr, Y_tr, Lambda, output_params, additional_arguments);
    
    cpu_time( ((t-1)*length(Lambda)*output_params.Kfolds + 1):(t*length(Lambda)*output_params.Kfolds) ) = cpu_time_tmp;
    
    % Caution: bias due to the CV !
    lambda_star = sqrt( (output_params.Kfolds-1)/output_params.Kfolds ) * lambda_star;
    
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Compute error 
    
    if strcmp( additional_arguments.model, 'islasso' ),
        
        J = eval_islasso_eta_trick( lambda_star );
        
    elseif strcmp( additional_arguments.model, 'slasso' ),
        
        J = eval_slasso_eta_trick( lambda_star );
     
    elseif strcmp( additional_arguments.model, 'lasso' ),
        
        J = eval_lasso( n_tr*lambda_star );
        
    elseif strcmp( additional_arguments.model, 'activeset_rectangle' ),
        
        J = eval_activeset_rectangle( lambda_star );
        
    elseif strcmp( additional_arguments.model, 'activeset' ),
        
        J = eval_activeset( lambda_star );
     
    elseif strcmp( additional_arguments.model, 'slasso_sq_SDPT3' ),
        
        J = eval_slasso_sq_SDPT3( lambda_star );
        
    else
        disp( '\n No model specified\n' );
        return
    end
         
    w_hat    = zeros(p,1);
    w_hat(J) = ( XtX_tr_over_n(J,J) )\( XtY_tr_over_n(J) );
    
    prediction_error(t) = (w - w_hat)' * XtX_test_over_n * (w - w_hat) / ( w' * XtX_test_over_n * w );
    
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Compute hull estimation error 
    
    G_Jc       = ~any(additional_arguments.G(J,:),1);
    JcHull_hat = any( mexAndLogic(~J, additional_arguments.G(:,G_Jc)), 2); 
    JHull_hat  = ~JcHull_hat;
    
    
    hull_estimation_error(t) = norm( JHull_hat - JHull, 2 );
    estimation_error(t)      = norm( w - w_hat ) / norm( w );
    
%==========================================================================

end

    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Nested functions
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function support = eval_islasso_eta_trick( lambda )
        
        support = true(p,1);
        
        for directions=1:size(additional_arguments.IDX,1),
            
            tmp_idx = additional_arguments.IDX(directions,1):additional_arguments.IDX(directions,2);
            
            tmp_w   = slasso_eta_trick(  X_tr, Y_tr, XtY_tr_over_n, XtX_tr_over_n, ...
                                        additional_arguments.G(:,tmp_idx), ... 
                                        additional_arguments.D(:,tmp_idx), lambda, [] );
 
            support = support & ( abs(tmp_w) > output_params.threshold );
            
        end
        
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function support = eval_slasso_eta_trick( lambda )

        tmp_w   = slasso_eta_trick(  X_tr, Y_tr, XtY_tr_over_n, XtX_tr_over_n, ...
                                    additional_arguments.G, ... 
                                    additional_arguments.D, lambda, [] );

        support = ( abs(tmp_w) > output_params.threshold );  
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function support = eval_lasso( lambda )
        
        lasso_params.mode       = 2;
        lasso_params.numThreads = 1;
        lasso_params.lambda     = lambda;
        
        tmp_w   = mexLasso( Y_tr, X_tr, lasso_params);

        support = ( abs(tmp_w) > output_params.threshold );  
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function support = eval_activeset_rectangle( lambda )
        
        activeset_params.maxCardJ       = output_params.maxCardJ;
        activeset_params.MaxIterations  = output_params.MaxIterations;
        
        tmp_w = activesetalgorithm_rectangle( X_tr, XtY_tr_over_n, ...
                                              additional_arguments.G, additional_arguments.D,...
                                              additional_arguments.IDX, lambda, height, width, activeset_params );

        support = ( abs(tmp_w) > output_params.threshold );  
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------  
    function support = eval_activeset( lambda )
        
        activeset_params.maxCardJ       = output_params.maxCardJ;
        activeset_params.MaxIterations  = output_params.MaxIterations;
        
        tmp_w   = activesetalgorithm( X_tr, XtY_tr_over_n, ...
                                      additional_arguments.G, additional_arguments.D, ...
                                      additional_arguments.IDX, lambda, activeset_params );

        support = ( abs(tmp_w) > output_params.threshold );  
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------  
    function support = eval_slasso_sq_SDPT3( lambda )
        
        SDPT3_params.squared_penalization_mode = true;
        
        tmp_w   = slasso_SDPT3(  XtX_tr_over_n, XtY_tr_over_n, ...
                                 additional_arguments.G, additional_arguments.D, lambda, SDPT3_params );

        support = ( abs(tmp_w) > output_params.threshold );  
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
end

%==========================================================================
%==========================================================================

function  [lambda_star, cpu_time] = get_lambda_by_CV( true_w, X_tr_, Y_tr_, Lambda, output_params, additional_arguments)

    non_empty_support       = true(length(Lambda),1);
    cpu_time                = zeros(output_params.Kfolds*length(Lambda), 1);
    CV_prediction_error     = zeros(length(Lambda),1);
    
    [ idx_train, idx_test ] = get_CV_idx( length( Y_tr_ ), output_params.Kfolds );
    
    for k=1:output_params.Kfolds,
        
        X_tr = X_tr_(idx_train{k},:);
        Y_tr = Y_tr_(idx_train{k},:);
        
        XtX_tr_over_n = X_tr' * X_tr / length( idx_train{k} );
        XtY_tr_over_n = X_tr' * Y_tr / length( idx_train{k} );
        
        X_test          = X_tr_(idx_test{k},:);
        XtX_test_over_n = X_test' * X_test / length( idx_test{k} );
        
        for j=1:length(Lambda),
            
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            
            tic
            
            if strcmp( additional_arguments.model, 'islasso' ),

                J = eval_islasso_eta_trick( Lambda(j) );

            elseif strcmp( additional_arguments.model, 'slasso' ),

                J = eval_slasso_eta_trick( Lambda(j) );

            elseif strcmp( additional_arguments.model, 'lasso' ),

                J = eval_lasso( length( idx_train{k} ) * Lambda(j) );% Caution : normalization for mexLasso !

            elseif strcmp( additional_arguments.model, 'activeset_rectangle' ),

                J = eval_activeset_rectangle( Lambda(j) );

            elseif strcmp( additional_arguments.model, 'activeset' ),

                J = eval_activeset( Lambda(j) );

            elseif strcmp( additional_arguments.model, 'slasso_sq_SDPT3' ),

                J = eval_slasso_sq_SDPT3( Lambda(j) ); 
                
            else
                
                disp( '\n No model specified\n' );
                return
                
            end
            
            cpu_time( (k-1)*length(Lambda)+j ) = toc;
            
            %--------------------------------------------------------------
            %--------------------------------------------------------------
            
            w_    = zeros(length(true_w),1);
            w_(J) = ( XtX_tr_over_n(J,J) )\( XtY_tr_over_n(J) );

            CV_prediction_error(j) = CV_prediction_error(j) + (w_-true_w)'*XtX_test_over_n*(w_-true_w) / (true_w'*XtX_test_over_n*true_w);
            
            % We check that for Lambda(j), there has never been an empty support
            non_empty_support(j) = non_empty_support(j) & any(J); 
            
        end
        
    end
   
    min_  = min( CV_prediction_error(non_empty_support) );
    
    % We take the sparsest model, i.e., the one obtained with the largest lambda
    
    lambda_star = Lambda( find( CV_prediction_error == min_, 1, 'last' ) );
    
%==========================================================================    
    
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Nested functions
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function support = eval_islasso_eta_trick( lambda )
        
        %------------------------------------------------------------------
        % warm start
        if ( j == 1 ), params = []; else params.w = w_; end;
        %------------------------------------------------------------------
        
        support = true( length(true_w), 1);
        
        for directions=1:size(additional_arguments.IDX,1),
            
            tmp_idx = additional_arguments.IDX(directions,1):additional_arguments.IDX(directions,2);
            
            tmp_w   = slasso_eta_trick(  X_tr, Y_tr, XtY_tr_over_n, XtX_tr_over_n, ...
                                        additional_arguments.G(:,tmp_idx), ... 
                                        additional_arguments.D(:,tmp_idx), lambda, params );
 
            support = support & ( abs(tmp_w) > output_params.threshold );
            
        end
        
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function support = eval_slasso_eta_trick( lambda )
        
        %------------------------------------------------------------------
        % warm start
        if ( j == 1 ), params = []; else params.w = w_; end;
        %------------------------------------------------------------------ 
        

        tmp_w   = slasso_eta_trick(  X_tr, Y_tr, XtY_tr_over_n, XtX_tr_over_n, ...
                                        additional_arguments.G, ... 
                                        additional_arguments.D, lambda, params );

        support = ( abs(tmp_w) > output_params.threshold );  
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    function support = eval_lasso( lambda )
        
        lasso_params.mode       = 2;
        lasso_params.numThreads = 1;
        lasso_params.lambda     = lambda;
        
        tmp_w   = mexLasso( Y_tr, X_tr, lasso_params);
        
        support = ( abs(tmp_w) > output_params.threshold );  
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------    
    function support = eval_activeset_rectangle( lambda )
        
        activeset_params.maxCardJ       = output_params.maxCardJ;
        activeset_params.MaxIterations  = output_params.MaxIterations;
        
        tmp_w = activesetalgorithm_rectangle( X_tr, XtY_tr_over_n, ...
                                              additional_arguments.G, additional_arguments.D,...
                                              additional_arguments.IDX, lambda, output_params.height, output_params.width, activeset_params );

        support = ( abs(tmp_w) > output_params.threshold );  
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------  
    function support = eval_activeset( lambda )
        
        activeset_params.maxCardJ       = output_params.maxCardJ;
        activeset_params.MaxIterations  = output_params.MaxIterations;
        
        tmp_w   = activesetalgorithm( X_tr, XtY_tr_over_n, ...
                                      additional_arguments.G, additional_arguments.D, ...
                                      additional_arguments.IDX, lambda, activeset_params );

        support = ( abs(tmp_w) > output_params.threshold );  
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------  
    function support = eval_slasso_sq_SDPT3( lambda )
        
        SDPT3_params.squared_penalization_mode = true;
        
        tmp_w   = slasso_SDPT3(  XtX_tr_over_n, XtY_tr_over_n, ...
                                 additional_arguments.G, additional_arguments.D, lambda, SDPT3_params );

        support = ( abs(tmp_w) > output_params.threshold );  
        
    end
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------

end

%==========================================================================
%==========================================================================

