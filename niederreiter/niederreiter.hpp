void calcc ( );
void calcv ( int px[], int b[], int v[], int v_max );
void golo ( double quasi[] );
int i4_characteristic ( int q );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
void inlo ( int dim, int base, int skip );
void niederreiter ( int dim_num, int base, int *seed, double r[] );
void niederreiter_generate ( int dim_num, int n, int base, int *seed, 
  double r[] );
void niederreiter_write ( int dim_num, int n, int base, int skip, double r[], 
  char *output_filename );
int *plymul ( int pa[], int pb[] );
double r8_epsilon ( );
void setfld ( int q );
void timestamp ( );

