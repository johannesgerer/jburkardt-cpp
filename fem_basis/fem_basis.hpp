double fem_basis_1d ( int i, int j, double x );
double fem_basis_2d ( int i, int j, int k, double x, double y );
double fem_basis_3d ( int i, int j, int k, int l, double x, double y, double z );
double fem_basis_md ( int m, int i[], double x[] );
double fem_basis_prism_triangle ( int i[], int j[], double xyz[] );
int i4vec_sum ( int n, int a[] );
double r8_fraction ( int i, int j );
double r8vec_sum ( int n, double a[] );
void timestamp ( );
