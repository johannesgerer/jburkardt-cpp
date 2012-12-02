# include <string>

namespace webbur
{
  void product_mixed_growth_weight ( int dim_num, int order_1d[], int order_nd, 
    int rule[], int np[], double p[], 
    void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ),
    double weight_nd[] );

  void sgmg_index ( int dim_num, int level_max, int rule[], 
    int point_num, int point_total_num, int sparse_unique_index[], int growth[],
    int sparse_order[], int sparse_index[] );

  void sgmg_point ( int dim_num, int level_max, int rule[], 
    int np[], double p[], 
    void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
    int point_num, int sparse_order[], int sparse_index[], int growth[],
    double sparse_point[] );

  int sgmg_size ( int dim_num, int level_max, int rule[], 
    int  np[], double p[], 
    void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
    double tol, int growth[] );

  int sgmg_size_total ( int dim_num, int level_max, int rule[],
    int growth[] );

  void sgmg_unique_index ( int dim_num, int level_max, int rule[], 
    int np[], double p[], 
    void ( *gw_compute_points[] ) ( int order, int np, double p[], double x[] ),
    double tol, int point_num, int point_total_num, int growth[],
    int sparse_unique_index[] );

  void sgmg_weight ( int dim_num, int level_max, int rule[], 
    int np[], double p[], 
    void ( *gw_compute_weights[] ) ( int order, int np, double p[], double w[] ),
    int point_num, int point_total_num, int sparse_unique_index[], int growth[],
    double sparse_weight[] );

  void sgmg_write ( int dim_num, int rule[], int np[], 
    double p[], int point_num, double sparse_weight[], double sparse_point[], 
    std::string file_name );
}
