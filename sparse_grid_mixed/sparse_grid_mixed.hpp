# include <string>

namespace webbur
{
  void sparse_grid_mixed_index ( int dim_num, int level_max, int rule[], 
    int point_num, int point_total_num, int sparse_unique_index[], 
    int sparse_order[], int sparse_index[] );

  void sparse_grid_mixed_point ( int dim_num, int level_max, int rule[], 
    double alpha[], double beta[], int point_num, int sparse_order[], 
    int sparse_index[], double sparse_point[] );

  int sparse_grid_mixed_size ( int dim_num, int level_max, int rule[], double alpha[],
    double beta[], double tol );

  int sparse_grid_mixed_size_total ( int dim_num, int level_max, int rule[] );

  void sparse_grid_mixed_unique_index ( int dim_num, int level_max, int rule[], 
    double alpha[], double beta[], double tol, int point_num, int point_total_num,
   int sparse_unique_index[] );

  void sparse_grid_mixed_weight ( int dim_num, int level_max, int rule[], 
    double alpha[], double beta[], int point_num, int point_total_num, 
    int sparse_unique_index[], double sparse_weight[] );

  void sparse_grid_mixed_write ( int dim_num, int rule[], double alpha[], 
    double beta[], int point_num, double sparse_weight[], double sparse_point[], 
    std::string file_name );
}
  
