int *catalan ( int n );
void catalan_values ( int &n_data, int &n, int &c );
void cbt_traverse ( int depth );
int graph_adj_edge_count ( int adj[], int nnode );
int graph_adj_is_node_connected ( int adj[], int nnode );
int graph_adj_is_tree ( int adj[], int nnode );
int *graph_arc_degree ( int nnode, int nedge, int inode[], int jnode[] );
int graph_arc_is_tree ( int nedge, int inode[], int jnode[], int nnode );
int graph_arc_node_count ( int nedge, int inode[], int jnode[] );
int graph_arc_node_max ( int nedge, int inode[], int jnode[] );
void graph_arc_print ( int nedge, int inode[], int jnode[], string title );
int *graph_arc_to_graph_adj ( int nedge, int inode[], int jnode[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_power ( int i, int j );
int i4_uniform_ab ( int a, int b, int &seed );
void i4mat_print ( int m, int n, int a[], string title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void i4vec_heap_d ( int n, int a[] );
int *i4vec_indicator_new ( int n );
int i4vec_max ( int n, int a[] );
void i4vec_print ( int n, int a[], string title );
void i4vec_sort_heap_a ( int n, int a[] );
int i4vec_sorted_unique_count ( int n, int a[] );
void pruefer_to_tree_arc ( int nnode, int iarray[], int inode[], int jnode[] );
void pruefer_to_tree_2 ( int nnode, int iarray[], int itree[] );
int *pruefer_to_tree_2_new ( int nnode, int iarray[] );
double r8_uniform_01 ( int &seed );
void timestamp ( );
void tree_arc_center ( int nnode, int inode[], int jnode[], int center[], 
  int &eccent, int &parity );
void tree_arc_diam ( int nnode, int inode[], int jnode[], int &diam, 
  int label[], int &n1, int &n2 );
void tree_arc_random ( int nnode, int &seed, int code[], int inode[], 
  int jnode[] );
int *tree_arc_to_pruefer ( int nnode, int inode[], int jnode[] );
int tree_enum ( int nnode );
void tree_parent_next ( int nnode, int code[], int itree[], int &more );
void tree_parent_to_arc ( int nnode, int parent[], int &nedge, int inode[], 
  int jnode[] );
int tree_rb_enum ( int n );
void tree_rb_lex_next ( int n, int a[], int &more );
int *tree_rb_to_parent ( int n, int a[] );
void tree_rb_yule ( int &n, int &seed, int a[] );
int *tree_rooted_code ( int nnode, int parent[] );
int tree_rooted_code_compare ( int nnode, int npart, int code1[], int code2[] );
void tree_rooted_depth ( int nnode, int parent[], int &depth, int depth_node[] );
int *tree_rooted_enum ( int nnode );
int *tree_rooted_random ( int nnode, int &seed );
void vec_next ( int n, int ibase, int iarray[], int &more );
void vec_random ( int n, int base, int &seed, int a[] );
int *vec_random_new ( int n, int base, int &seed );
