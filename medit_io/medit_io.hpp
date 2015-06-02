char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
void cyl248_data ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[],
  int hexahedron_label[] );
void cyl248_size ( int *dim, int *vertices, int *edges, int *triangles,
  int *quadrilaterals, int *tetrahedrons, int *hexahedrons );
void hexahexa_2x2x2_data ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[], int hexahedron_vertex[],
  int hexahedron_label[] );
void hexahexa_2x2x2_size ( int *dim, int *vertices, int *edges, int *triangles,
  int *quadrilaterals, int *tetrahedrons, int *hexahedrons );
void i4vec_copy ( int n, int a1[], int a2[] );
void i4vec_zero ( int n, int a[] );
void mesh_data_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void mesh_data_read ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void mesh_size_print ( int dim, int vertices, int edges, int triangles,
  int quadrilaterals, int tetrahedrons, int hexahedrons );
void mesh_size_read ( string filename, int *dim, int *vertices, int *edges,
  int *triangles, int *quadrilaterals, int *tetrahedrons, int *hexahedrons );
void mesh_write ( string filename, int dim, int vertices, int edges,
  int triangles, int quadrilaterals, int tetrahedrons, int hexahedrons,
  double vertex_coordinate[], int vertex_label[], int edge_vertex[],
  int edge_label[], int triangle_vertex[], int triangle_label[],
  int quadrilateral_vertex[], int quadrilateral_label[],
  int tetrahedron_vertex[], int tetrahedron_label[],
  int hexahedron_vertex[], int hexahedron_label[] );
void r8vec_copy ( int n, double a1[], double a2[] );
void r8vec_zero ( int n, double a[] );
bool s_begin ( string s1, string s2 );
bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
string s_newline_to_null ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
void timestamp ( );
