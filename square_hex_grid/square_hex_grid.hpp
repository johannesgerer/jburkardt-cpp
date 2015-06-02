void box_print_2d ( double box[] );
char digit_to_ch ( int i );
void hex_grid_01_approximate_h ( double h_goal, int *nodes_per_layer, 
  double *h );
void hex_grid_01_approximate_n ( int n_goal, int *nodes_per_layer, int *n );
void hex_grid_01_h ( int nodes_per_layer, double *hx, double *hy );
int hex_grid_01_layers ( int nodes_per_layer );
int hex_grid_01_n ( int nodes_per_layer );
double *hex_grid_01_points ( int nodes_per_layer, int layers, int n );
void hex_grid_approximate_h ( double box[], double h_goal, 
  int *nodes_per_layer, double *h );
void hex_grid_approximate_n ( double box[], int n_goal, int *nodes_per_layer, 
  int *n );
void hex_grid_h ( int nodes_per_layer, double box[], double *hx, double *hy );
int hex_grid_layers ( int nodes_per_layer, double box[] );
int hex_grid_n ( int nodes_per_layer, double box[] );
double *hex_grid_points ( int nodes_per_layer, int layers, int n, 
  double box[] );
void hex_grid_write ( int n, int nodes_per_layer, int layers, double hx, 
  double hy, double box[], double r[], char *file_out_name );
int i4_log_10 ( int i );
char *i4_to_s ( int i );
double r8_epsilon ( void );
void timestamp ( void );
char *timestring ( void );
