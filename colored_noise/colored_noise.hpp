double *f_alpha ( int n, double q_d, double alpha, int *seed );
double r8_normal_01 ( int *seed );
string r8_to_string ( double r8, string format );
double r8_uniform_01 ( int *seed );
void r8vec_print ( int n, double a[], string title );
void r8vec_print_part ( int n, double a[], int max_print, string title );
double *r8vec_sftb ( int n, double azero, double a[], double b[] );
void r8vec_sftf ( int n, double r[], double *azero, double a[], double b[] );
void timestamp ( );
