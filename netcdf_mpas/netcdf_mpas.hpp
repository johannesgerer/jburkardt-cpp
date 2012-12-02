void netcdf_mpas_read_cellsonedge ( string filename, int ncells, 
  int cellsonedge[] );
void netcdf_mpas_read_cellsonvertex ( string filename, int nvertices,
  int cellsonvertex[] );
int netcdf_mpas_read_maxedges ( string filename );
int netcdf_mpas_read_ncells ( string filename );
int netcdf_mpas_read_nedges ( string filename );
void netcdf_mpas_read_nedgesoncell ( string filename, int ncells, 
  int nedgesoncell[] );
int netcdf_mpas_read_nvertices ( string filename );
void netcdf_mpas_read_verticesoncell ( string filename, int maxedges,
  int ncells, int verticesoncell[] );
void netcdf_mpas_read_xyzcell ( string filename, int ncells, double xcell[], 
  double ycell[], double zcell[] );
void netcdf_mpas_read_xyzvertex ( string filename, int nvertices,
  double xvertex[], double yvertex[], double zvertex[] );
void netcdf_mpas_report ( string filename );
void timestamp ( );
