# include <string>

namespace webbur
{
  void product_mixed_growth_weight 
  ( 
    int dim_num, 
    int order_1d[], 
    int order_nd, 
    void ( *gw_compute_weights[] ) ( int order, int dim, double w[] ), 
    double weight_nd[] 
  );

  void sgmg_index
  ( 
    int dim_num,
    int level_max, 
    int point_num, 
    int point_total_num, 
    int sparse_unique_index[], 
    int growth,
    int ( *gw_compute_order[] ) ( int level, int growth ),
    int sparse_order[], 
    int sparse_index[] 
  );

  void sgmg_point
  ( 
    int dim_num, 
    int level_max,
    void ( *gw_compute_points[] ) ( int order, int dim, double x[] ),
    int point_num,
    int sparse_order[],
    int sparse_index[],
    int growth,
    int ( *gw_compute_order[] ) ( int level, int growth ),
    double sparse_point[]
  );

  int sgmg_size 
  ( 
    int dim_num, 
    int level_max, 
    void ( *gw_compute_points[] ) ( int order, int dim, double x[] ),
    double tol, 
    int growth,
    int ( *gw_compute_order[] ) ( int level, int growth )  
  );

  int sgmg_size_total
  ( 
    int dim_num,
    int level_max,
    int growth,
    int ( *gw_compute_order[] ) ( int level, int growth ) 
  );

  void sgmg_unique_index
  ( 
    int dim_num, 
    int level_max,
    void ( *gw_compute_points[] ) ( int order, int dim, double x[] ),
    double tol, 
    int point_num,
    int point_total_num,
    int growth, 
    int ( *gw_compute_order[] ) ( int level, int growth ),
    int sparse_unique_index[] 
  );

  void sgmg_weight
  ( 
    int dim_num, 
    int level_max,
    void ( *gw_compute_weights[] ) ( int order, int dim, double w[] ),
    int point_num,
    int point_total_num,
    int sparse_unique_index[],
    int growth,
    int ( *gw_compute_order[] ) ( int level, int growth ),
    double sparse_weight[] 
  );

}
