#include <string>
#include <netcdfcpp.h>

int netcdf_mpas_read_dim ( string filename, string dim_name );
int netcdf_mpas_read_num_vars(string filename);
void netcdf_mpas_read_xyzcell ( string filename, int ncells, double xcell[], double ycell[], double zcell[] );
void netcdf_mpas_read_xyzvertex ( string filename, int nvertices, double xvertex[], double yvertex[], double zvertex[] );
void netcdf_mpas_read_xyzedge ( string filename, int nedges, double xedge[], double yedge[], double zedge[] );
void netcdf_mpas_read_verticesoncell ( string filename, int maxedges, int ncells, int verticesOnCell[] );
void netcdf_mpas_read_verticesonedge ( string filename, int nedges, int verticesonedge[] );
void netcdf_mpas_read_cellsonvertex ( string filename, int nvertices, int cellsonvertex[] );
void netcdf_mpas_read_cellsonedge ( string filename, int nedges, int cellsonedge[] );
void netcdf_mpas_read_nedgesoncell ( string filename, int ncells, int nedgesoncell[] );
int netcdf_mpas_list_ncell_fields(string filename);
int netcdf_mpas_list_nvertex_fields(string filename);
int netcdf_mpas_list_nedge_fields(string filename);
int netcdf_mpas_field_num_dims(string filename, int id);
int netcdf_mpas_field_num_items(string filename, int id);
void netcdf_mpas_print_field_info(string filename, int id);
void netcdf_mpas_read_field(string filename, int id, double values[], int cur_time, int cur_level);
void netcdf_mpas_read_full_field(string filename, int id, double values[]);
