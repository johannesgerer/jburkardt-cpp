int pyramid_grid_size ( int n );
double *pyramid_unit_grid ( int n, int ng );
void pyramid_unit_grid_plot ( int n, int ng, double pg[], string header );
void pyramid_unit_vertices ( double v1[], double v2[], double v3[], 
  double v4[], double v5[] );
void r8_print ( double r, string title );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void timestamp ( );

