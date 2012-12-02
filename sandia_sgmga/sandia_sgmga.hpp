# include <string>

namespace webbur
{
  double *sandia_sgmga_aniso_balance 
  ( 
    double alpha_max,
    int dim_num, 
    double level_weight[] 
  );

  void sandia_sgmga_aniso_normalize 
  ( 
    int option,
    int dim_num, 
    double level_weight[]
  );

  void sandia_sgmga_importance_to_aniso 
  ( 
    int dim_num, 
    double importance[], 
    double level_weight[]
  );

  void sandia_sgmga_index 
  (
    int dim_num, 
    double level_weight[], 
    int level_max, 
    int point_num, 
    int point_total_num, 
    int sparse_unique_index[],
    int growth,
    int ( *gw_compute_order[] ) ( int level, int growth ),
    int sparse_order[], 
    int sparse_index[]
  );

  void sandia_sgmga_point 
  ( 
    int dim_num, 
    double level_weight[],
    int level_max, 
    void ( *gw_compute_points[] ) ( int order, int dim, double x[] ),
    int point_num,
    int sparse_order[],
    int sparse_index[], 
    int growth, 
    int ( *gw_compute_order[] ) ( int level, int growth ),
    double sparse_point[]
  );

  void sandia_sgmga_product_weight 
  (
    int dim_num, 
    int order_1d[], 
    int order_nd,  
    void ( *gw_compute_weights[] ) ( int order, int dim, double w[] ),
    double weight_nd[] 
  );

  int sandia_sgmga_size 
  (
    int dim_num,
    double level_weight[],
    int level_max, 
    void ( *gw_compute_points[] ) ( int order, int dim, double x[] ),
    double tol,
    int growth,
    int ( *gw_compute_order[] ) ( int level, int growth ) 
  );

  int sandia_sgmga_size_total 
  ( 
    int dim_num,
    double level_weight[],
    int level_max, 
    int growth,
    int ( *gw_compute_order[] ) ( int level, int growth ) 
  );

  void sandia_sgmga_unique_index 
  ( 
    int dim_num,
    double level_weight[], 
    int level_max,
    void ( *gw_compute_points[] ) ( int order, int dim, double x[] ),
    double tol,
    int point_num,
    int point_total_num,
    int growth, 
    int ( *gw_compute_order[] ) ( int level, int growth ),
    int sparse_unique_index[] 
  );

  void sandia_sgmga_vcn 
  (
    int n,
    double level_weight[],
    int x[], 
    double q_min,
    double q_max, 
    bool *more 
  );

  double sandia_sgmga_vcn_coef 
  (
    int n, 
    double level_weight[],
    int x[], 
    double q_max 
  );

  double sandia_sgmga_vcn_coef_naive 
  (
    int n,
    double level_weight[],
    int x_max[],
    int x[],
    double q_min,
    double q_max 
  );

  void sandia_sgmga_vcn_naive 
  (
    int n,
    double level_weight[],
    int x_max[],
    int x[], 
    double q_min,
    double q_max,
    bool *more 
  );

  void sandia_sgmga_vcn_ordered 
  (
    int dim_num,
    double level_weight[],
    int x_max[],
    int x[],
    double q_min,
    double q_max,
    bool *more 
  );

  void sandia_sgmga_vcn_ordered_naive 
  (
    int dim_num,
    double level_weight[],
    int x_max[],
    int x[],
    double q_min,
    double q_max,
    bool *more 
  );

  void sandia_sgmga_weight 
  (
    int dim_num,
    double level_weight[],
    int level_max, 
    void ( *gw_compute_weights[] ) ( int order, int dim, double w[] ),
    int point_num, 
    int point_total_num,
    int sparse_unique_index[], 
    int growth,
    int ( *gw_compute_order[] ) ( int level, int growth ),
    double sparse_weight[] );
}
