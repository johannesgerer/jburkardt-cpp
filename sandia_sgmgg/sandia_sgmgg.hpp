namespace webbur
{
  void sandia_sgmgg_coef_inc2 ( int m, int n1, int s1[], int c1[],
    int s2[], int c3[] );
  void sandia_sgmgg_coef_naive ( int dim_num, int point_num, int sparse_index[],
    int coef[] );
  void sandia_sgmgg_neighbor_naive ( int dim_num, int point_num, int sparse_index[],
    int neighbor[] );
}
