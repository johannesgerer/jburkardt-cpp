# include <fstream>
# include <string>

namespace webbur
{
  void ccn_points ( int n, int dim, double x[] );
  void ccn_weights ( int n, int dim, double w[] );

  void clenshaw_curtis_points ( int n, int dim, double x[] );
  void clenshaw_curtis_weights ( int n, int dim, double w[] );

  void fejer2_points ( int n, int dim, double x[] );
  void fejer2_weights ( int n, int dim, double w[] );

  void gen_hermite_points ( int n, int dim, double x[] );
  void gen_hermite_weights ( int n, int dim, double w[] );

  void gen_laguerre_points ( int n, int dim, double x[] );
  void gen_laguerre_weights ( int n, int dim, double w[] );

  void hcc_points ( int n, int dim, double x[] );
  void hcc_weights ( int n, int dim, double w[] );

  void hce_points ( int n, int dim, double x[] );
  void hce_weights ( int n, int dim, double w[] );

  void hermite_genz_keister_points ( int n, int dim, double x[] );
  void hermite_genz_keister_weights ( int n, int dim, double w[] );

  void hermite_points ( int n, int dim, double x[] );
  void hermite_weights ( int n, int dim, double w[] );

  void jacobi_points ( int n, int dim, double x[] );
  void jacobi_weights ( int n, int dim, double w[] );

  void laguerre_points ( int n, int dim, double x[] );
  void laguerre_weights ( int n, int dim, double w[] );

  void legendre_points ( int n, int dim, double x[] );
  void legendre_weights ( int n, int dim, double w[] );

  void ncc_points ( int n, int dim, double x[] );
  void ncc_weights ( int n, int dim, double w[] );

  void nco_points ( int n, int dim, double x[] );
  void nco_weights ( int n, int dim, double w[] );

  double parameter ( int dim, int offset );

  void patterson_points ( int n, int dim, double x[] );
  void patterson_weights ( int n, int dim, double w[] );
}

