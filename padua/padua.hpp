void filename_inc ( string *filename );
string i4_to_string ( int i4 );
int padua_order ( int l );
void padua_plot ( int l, string filename );
void padua_points_set ( int l, double x[], double y[] );
double *padua_points ( int l );
double *padua_weights_set ( int l );
double *padua_weights ( int l );
double r8_max ( double x, double y );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void r8vec_reverse ( int n, double a[] );
void timestamp ( );

