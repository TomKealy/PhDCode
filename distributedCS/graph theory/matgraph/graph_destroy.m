function graph_destroy
% destroy the GRAPH_MAGIC system

disp('WARNING: This will delete all graphs and invalidate all graph handles')
yn = input('Are you sure you want to proceed? 1 = YES, 0 = NO --> ');

if (yn == 1)
    clear global GRAPH_MAGIC
    disp('All graphs destroyed. Start over with graph_init');
else
    disp('Destruction canceled');
end
    