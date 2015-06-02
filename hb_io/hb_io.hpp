bool ch_eqi ( char c1, char c2 );
bool ch_is_digit ( char c );
bool ch_is_format_code ( char c );
int ch_to_digit ( char c );
void hb_exact_read ( ifstream &input, int nrow, int nrhs, int rhscrd, 
  char *rhsfmt, char *rhstyp, double exact[] );
void hb_exact_write ( ofstream &output, int nrow, int nrhs, int rhscrd, 
  char *rhsfmt, char *rhstyp, double exact[] );
void hb_file_read ( ifstream &input, char **title, char **key, int *totcrd, 
  int *ptrcrd, int *indcrd, int *valcrd, int *rhscrd, char **mxtype, int *nrow, 
  int *ncol, int *nnzero, int *neltvl, char **ptrfmt, char **indfmt, char **valfmt, 
  char **rhsfmt, char **rhstyp, int *nrhs, int *nrhsix, int **colptr, 
  int **rowind, double **values, double **rhsval, int **rhsptr, int **rhsind,  
  double **rhsvec, double **guess, double **exact );
void hb_file_write ( ofstream &output, char *title, char *key, int totcrd, 
  int ptrcrd, int indcrd, int valcrd, int rhscrd, char *mxtype, int nrow, 
  int ncol, int nnzero, int neltvl, char *ptrfmt, char *indfmt, char *valfmt, 
  char *rhsfmt, char *rhstyp, int nrhs, int nrhsix, int colptr[],
  int rowind[], double values[], double rhsval[], int rhsptr[], int rhsind[], 
  double rhsvec[], double guess[], double exact[] );
void hb_guess_read ( ifstream &input, int nrow, int nrhs, int rhscrd, 
  char *rhsfmt, char *rhstyp, double guess[] );
void hb_guess_write ( ofstream &output, int nrow, int nrhs, int rhscrd, 
  char *rhsfmt, char *rhstyp, double guess[] );
void hb_header_print ( char *title, char *key, int totcrd, int ptrcrd, 
  int indcrd, int valcrd, int rhscrd, char *mxtype, int nrow, int ncol, 
  int nnzero, int neltvl, char *ptrfmt, char *indfmt, char *valfmt, 
  char *rhsfmt, char *rhstyp, int nrhs, int nrhsix );
void hb_header_read ( ifstream &input, char **title, char **key, int *totcrd, 
  int *ptrcrd, int *indcrd, int *valcrd, int *rhscrd, char **mxtype, int *nrow, 
  int *ncol, int *nnzero, int *neltvl, char **ptrfmt, char **indfmt, 
  char **valfmt, char **rhsfmt, char **rhstyp, int *nrhs, int *nrhsix );
void hb_header_write ( ofstream &output, char *title, char *key, int totcrd, 
  int ptrcrd, int indcrd, int valcrd, int rhscrd, char *mxtype, int nrow, 
  int ncol, int nnzero, int neltvl, char *ptrfmt, char *indfmt, char *valfmt, 
  char *rhsfmt, char *rhstyp, int nrhs, int nrhsix );
double *hb_matvec_a_mem ( int nrow, int ncol, int nnzero, int nrhs, 
  int colptr[], int rowind[], double values[], double exact[] );
void hb_rhs_read ( ifstream &input, int nrow, int nnzero, int nrhs, int nrhsix, 
  int rhscrd, char *ptrfmt, char *indfmt, char *rhsfmt, char *mxtype, 
  char *rhstyp, double rhsval[], int rhsind[], int rhsptr[], double rhsvec[] );
void hb_rhs_write ( ofstream &output, int nrow, int nnzero, int nrhs, int nrhsix, 
  int rhscrd, char *ptrfmt, char *indfmt, char *rhsfmt, char *mxtype, 
  char *rhstyp, double rhsval[], int rhsind[], int rhsptr[], double rhsvec[] );
void hb_structure_print ( int ncol, char *mxtype, int nnzero, int neltvl, 
  int colptr[], int rowind[] );
void hb_structure_read ( ifstream &input, int ncol, char *mxtype, int nnzero, 
  int neltvl, int ptrcrd, char *ptrfmt, int indcrd, char *indfmt, 
  int colptr[], int rowind[] );
void hb_structure_write ( ofstream &output, int ncol, char *mxtype, 
  int nnzero, int neltvl, char *ptrfmt, char *indfmt, int colptr[], 
  int rowind[] );
int *hb_ua_colind ( int ncol, int colptr[], int nnzero );
void hb_values_print ( int ncol, int colptr[], char *mxtype, int nnzero, 
  int neltvl, double values[] );
void hb_values_read ( ifstream &input, int valcrd, char *mxtype, int nnzero,
  int neltvl, char *valfmt, double values[] );
void hb_values_write ( ofstream &output, int valcrd, char *mxtype, 
  int nnzero, int neltvl, char *valfmt, double values[] );
double *hb_vecmat_a_mem ( int nrow, int ncol, int nnzero, int nrhs, 
  int colptr[], int rowind[], double values[], double exact[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4vec_print ( int n, int a[], string title );
void i4vec_print_part ( int n, int a[], int max_print, string title );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
void r8vec_print ( int n, double a[], string title );
void r8vec_print_part ( int n, double a[], int max_print, string title );
int s_len_trim ( char *s );
char *s_substring ( char *s, int a, int b );
void s_to_format ( char *s, int *r, char *code, int *w, int *m );
void s_trim ( char *s );
void timestamp ( );
