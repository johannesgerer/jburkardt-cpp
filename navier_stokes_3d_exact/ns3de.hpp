double r8vec_amax ( int n, double a[] );
double r8vec_amin ( int n, double a[] );
double r8vec_max ( int n, double r8vec[] );
double r8vec_min ( int n, double r8vec[] );
double *r8vec_uniform_ab_new ( int n, double a, double b, int &seed );
void resid_ethier ( double a, double d, int n, double x[], double y[], 
  double z[], double t, double ur[], double vr[], double wr[], double pr[] );
void timestamp ( );
void uvwp_ethier ( double a, double d, int n, double x[], double y[], 
  double z[], double t, double u[], double v[], double w[], double p[] );
