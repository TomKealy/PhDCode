function get_SDPT3_pi4_speed_comparison( p )

n_tr   = 3500;
n_runs = 5;

switch p
    case 100
        
        hull_pattern = [3:4, 12:15, 21:26, 31:36, 42:45, 53:54];
        height       = 10;
        width        = 10;
        
    case 225
        
        hull_pattern = [3:4, 17:20, 31:36, 46:51, 62:65, 78:79];
        height       = 15;
        width        = 15;
        
    case 400
        
        hull_pattern = [ 3:4, 22:25, 41:46, 61:66, 82:85, 103:104 ];
        height       = 20;
        width        = 20;   
        
    case 900
        
        hull_pattern = [ 3:4, 32:35, 61:66, 91:96, 122:125, 153:154 ];
        height       = 30;
        width        = 30;   
         
    case 1600
        
        hull_pattern = [ 3:4, 42:45, 81:86, 121:126, 162:165, 203:204 ];
        height       = 40;
        width        = 40;
        
    case 2500
        
        hull_pattern = [ 3:4, 52:55, 101:106, 151:156, 202:205, 253:254 ];
        height       = 50;
        width        = 50;
end


 params.slope     = 0;
 [ G1, D1 ] = get_groups( height, width, params );

 params.slope     = 1;
 [ G2, D2 ] = get_groups( height, width, params );

 G   = [G1, G2];
 D   = [D1, D2];
 

 Lambda = ( 2.^(-5:0.2:5) )/sqrt(n_tr);
 

 additional_arguments.G     = G;
 additional_arguments.D     = D;
 additional_arguments.model = 'slasso_sq_SDPT3'; 
 
[ prediction_error, estimation_error, hull_estimation_error, cpu_time, output_params ] = get_regpaths( height, width, hull_pattern, Lambda, n_tr, n_runs, additional_arguments );
 
 path_to_data = '~/OUTPUT_MATLAB/';
 
 save([path_to_data 'SDPT3_pi4_speed_comparison_p_' int2str(p)],'prediction_error', 'estimation_error', 'hull_estimation_error', 'cpu_time', 'output_params');

  
  
 
 
 
