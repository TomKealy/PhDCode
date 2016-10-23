function  get_simulation_slasso_sequence_exp( n_tr, convex_rate ) 

 n_runs = 50;

 height = 400;
 width  = 1;
 
  switch convex_rate
     
    case 25
        hull_pattern = [ 1:3, 22:24 ];
 
    case 50
        hull_pattern = [ 1:6, 19:24];
        
    case 75
        hull_pattern = [ 1:9, 16:24];
        
    case 100
        hull_pattern = 1:24;
 end
 
 [G,D] = get_groups( height, width, [] );

 Lambda = (2.^(-3:0.2:5) )/sqrt(n_tr);
 
 
 additional_arguments.G     = G;
 additional_arguments.D     = D;
 additional_arguments.model = 'slasso'; 
 
[ prediction_error, estimation_error, hull_estimation_error, cpu_time, output_params ] = get_regpaths( height, width, hull_pattern, Lambda, n_tr, n_runs, additional_arguments );

 
  path_to_data = '~/OUTPUT_MATLAB/';
 
 save([path_to_data 'slasso_sequence_exp_n_' int2str(n_tr) '_convex_' int2str(convex_rate)], 'prediction_error', 'estimation_error', 'hull_estimation_error', 'cpu_time', 'output_params');

end
