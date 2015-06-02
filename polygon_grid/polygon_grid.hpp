int polygon_grid_count ( int n, int nv );
void polygon_grid_display ( int n, int nv, double v[], int ng, double xg[], 
  string prefix );
double *polygon_grid_points ( int n, int nv, double v[], int ng );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );

