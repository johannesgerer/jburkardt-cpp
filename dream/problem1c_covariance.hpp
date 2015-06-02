//****************************************************************************80

class Covariance

//****************************************************************************80
{
  public:
    Covariance ( int par_num );
    ~Covariance ( );
    double *array_get ( );
    void array_set ( );
    double det_get ( );
    void det_set ( );
    double *factor_get ( );
    void factor_set ( );
    double *inv_get ( );
    void inv_set ( );
    double *mean_get ( );
    void mean_set ( );
    void print ( string title );

  private:
    int order;
    double *array;
    double det;
    double *factor;
    double *inv;
    double *mean;
};

double *r8mat_copy_new ( int m, int n, double a1[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
double *r8vec_copy_new ( int n, double a1[] );
