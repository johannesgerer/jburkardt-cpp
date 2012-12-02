//  blend.h

void blend_i_0d1   ( double x[], int m );

void blend_ij_0d1  ( double x[], int m1, int m2 );
void blend_ij_1d1  ( double x[], int m1, int m2 );

void blend_ijk_0d1 ( double x[], int m1, int m2, int m3 );
void blend_ijk_1d1 ( double x[], int m1, int m2, int m3 );
void blend_ijk_2d1 ( double x[], int m1, int m2, int m3 );

void blend_0d1     ( double r, double x0, double x1, double *x );

void blend_1d1     ( double r, double s, double x00, double x01, double x10, 
                     double x11, double xr0, double xr1, double x0s, 
                     double x1s, double *x );

void blend_2d1     ( double r, double s, double t, 
                     double x000, double x001, double x010, double x011, 
                     double x100, double x101, double x110, double x111, 
                     double xr00, double xr01, double xr10, double xr11, 
                     double x0s0, double x0s1, double x1s0, double x1s1, 
                     double x00t, double x01t, double x10t, double x11t, 
                     double x0st, double x1st, double xr0t, 
                     double xr1t, double xrs0, double xrs1, 
                     double *x );

void blend_r_0dn   ( double r, double x[], int n, 
                     void ( *bound_r ) ( double r, int i, double *xi ) );

void blend_rs_0dn  ( double r, double s, double x[], int n, 
                     void ( *bound_rs ) ( double r, double s, int i, double *xi ) );

void blend_rs_1dn  ( double r, double s, double x[], int n, 
                     void ( *bound_rs ) ( double r, double s, int i, double *xi ) );

void blend_rst_0dn ( double r, double s, double t, double x[], int n, 
                     void ( *bound_rst ) ( double r, double s, double t, int i, 
                     double *xi ) );

void blend_rst_1dn ( double r, double s, double t, double x[], int n, 
                     void ( *bound_rst ) ( double r, double s, double t, int i, 
                     double *xi ) );

void blend_rst_2dn ( double r, double s, double t, double x[], int n, 
                     void ( *bound_rst ) ( double r, double s, double t, int i, 
                     double *xi ) );

int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void r8block_print ( int l, int m, int n, double a[], char *title );
void r8mat_print ( int m, int n, double a[], char *title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
int s_len_trim ( char *s );
void timestamp ( void );
