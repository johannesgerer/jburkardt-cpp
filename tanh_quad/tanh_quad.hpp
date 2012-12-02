int i4_power ( int i, int j );
void midpoint_rule ( int n, double x[], double w[] );
double r8_abs ( double x );
double r8_epsilon ( );
double r8vec_dot ( int n, double a1[], double a2[] );
double r8vec_sum ( int n, double a[] );
void rule_adjust ( double a, double b, double c, double d, int order, 
  double x[], double w[] );
int tanh_h_to_n ( double h, double tol );
double tanh_m_to_h ( int m );
double tanh_n_to_h ( int n  );
void tanh_rule ( int n, double h, double x[], double w[] );
int tanh_sinh_h_to_n ( double h, double tol );
void tanh_sinh_rule ( int n, double h, double x[], double w[] );
void timestamp ( void );
void trap_rule ( int n, double x[], double w[] );
