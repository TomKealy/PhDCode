function  get_simulation_islasso_grid_exp_pi4( n_tr, convex_rate )
%GET_SIMULATION_LASSO_SEQUENCES 

 n_runs = 50;

 height = 20;
 width  = 20;
 
 switch convex_rate
     
    case 25
        hull_pattern = [ 3,4,41,46,61,66,103,104 ];
 
    case 50
        hull_pattern = [ 3,4,41,46,61,66,103,104,22,25,82,85];
        
    case 75
        hull_pattern = setdiff( [ 3:4, 22:25, 41:46, 61:66, 82:85, 103:104 ], [ 43, 44, 63, 64 ] );
        
    case 100
        hull_pattern = [ 3:4, 22:25, 41:46, 61:66, 82:85, 103:104 ];
 end
 
 params.intersection = true;
 params.slope     = 0;
 [ G1, D1, IDX1 ] = get_groups( height, width, params );

 params.slope     = 1;
 [ G2, D2, IDX2 ] = get_groups( height, width, params );

 G   = [G1, G2];
 D   = [D1, D2];
 IDX = [IDX1; IDX2+max(IDX1(:))];

 Lambda = (2.^(-5:0.2:5) )/sqrt(n_tr);
 
 
 additional_arguments.G     = G;
 additional_arguments.D     = D;
 additional_arguments.IDX   = IDX;
 additional_arguments.model = 'islasso'; 
 
[ prediction_error, estimation_error, hull_estimation_error, cpu_time, output_params ] = get_regpaths( height, width, hull_pattern, Lambda, n_tr, n_runs, additional_arguments );
 
 path_to_data = '~/OUTPUT_MATLAB/';
 
 save([path_to_data 'islasso_grid_exp_pi4_n_' int2str(n_tr) '_convex_' int2str(convex_rate)],'prediction_error', 'estimation_error', 'hull_estimation_error', 'cpu_time', 'output_params');


end
