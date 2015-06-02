void kjacopols ( double x, double a, double b, int n, double pols[] );
void kjacopols2 ( double x, double a, double b, int n, double pols[], 
  double ders[] );
double *klegeypols ( double x, double y, int n );
void klegeypols3 ( double x, double y, int n, double pols[], double dersx[], 
  double dersy[] );
double *ortho2eva ( int mmax, double z[] );
double *ortho2eva0 ( int mmax, double z[] );
void quaecopy ( int kk, double xs[], double ys[], double ws[], double z[], double w[] );
void quaecopy2 ( double xs[], double ys[], double ws[], double xnew[], 
  double ynew[], double w[], int kk );
int quaeinside ( int iitype, double xsout, double ysout );
void quaenodes ( int nptsout, double xsout[], double ysout[], double wsout[], 
  int nptsoutout, double xs2[], double ys2[], double ws2[] );
void quaenodes2 ( int nptsout, double xsout[], double ysout[], double wsout[], 
  int nptsoutout, double xs2[], double ys2[], double ws2[] );
void quaequad ( int itype, int mmax, double zs[], double whts[], int numnodes );
void quaequad0 ( int mmax, int kk, double xnew[], double ynew[], double w[] );
void quaerotate ( double xin, double yin, double &xout, double &yout );
void r8vec_copy ( int n, double a1[], double a2[] );
double r8vec_sum ( int n, double a[] );
double *r8vec_uniform_01_new ( int n, int &seed );
double *ref_to_koorn ( double r[] );
double *ref_to_triangle ( double tvert1[], double tvert2[], double tvert3[], 
  double r[] );
int rule_compressed_size ( int mmax );
int rule_full_size ( int mmax );
void rule01 ( double x[], double y[], double w[] );
void rule02 ( double x[], double y[], double w[] );
void rule03 ( double x[], double y[], double w[] );
void rule04 ( double x[], double y[], double w[] );
void rule05 ( double x[], double y[], double w[] );
void rule06 ( double x[], double y[], double w[] );
void rule07 ( double x[], double y[], double w[] );
void rule08 ( double x[], double y[], double w[] );
void rule09 ( double x[], double y[], double w[] );
void rule10 ( double x[], double y[], double w[] );
void rule11 ( double x[], double y[], double w[] );
void rule12 ( double x[], double y[], double w[] );
void rule13 ( double x[], double y[], double w[] );
void rule14 ( double x[], double y[], double w[] );
void rule15 ( double x[], double y[], double w[] );
void rule16 ( double x[], double y[], double w[] );
void rule17 ( double x[], double y[], double w[] );
void rule18 ( double x[], double y[], double w[] );
void rule19 ( double x[], double y[], double w[] );
void rule20 ( double x[], double y[], double w[] );
void rule21 ( double x[], double y[], double w[] );
void rule22 ( double x[], double y[], double w[] );
void rule23 ( double x[], double y[], double w[] );
void rule24 ( double x[], double y[], double w[] );
void rule25 ( double x[], double y[], double w[] );
void rule26 ( double x[], double y[], double w[] );
void rule27 ( double x[], double y[], double w[] );
void rule28 ( double x[], double y[], double w[] );
void rule29 ( double x[], double y[], double w[] );
void rule30 ( double x[], double y[], double w[] );
void rule31 ( double x[], double y[], double w[] );
void rule32 ( double x[], double y[], double w[] );
void rule33 ( double x[], double y[], double w[] );
void rule34 ( double x[], double y[], double w[] );
void rule35 ( double x[], double y[], double w[] );
void rule36 ( double x[], double y[], double w[] );
void rule37 ( double x[], double y[], double w[] );
void rule38 ( double x[], double y[], double w[] );
void rule39 ( double x[], double y[], double w[] );
void rule40 ( double x[], double y[], double w[] );
void rule41 ( double x[], double y[], double w[] );
void rule42 ( double x[], double y[], double w[] );
void rule43 ( double x[], double y[], double w[] );
void rule44 ( double x[], double y[], double w[] );
void rule45 ( double x[], double y[], double w[] );
void rule46 ( double x[], double y[], double w[] );
void rule47 ( double x[], double y[], double w[] );
void rule48 ( double x[], double y[], double w[] );
void rule49 ( double x[], double y[], double w[] );
void rule50 ( double x[], double y[], double w[] );
double *simplex_to_triangle ( double tvert1[], double tvert2[], double tvert3[], 
  double s[] );
void timestamp ( );
double triangle_area ( double vert1[], double vert2[], double vert3[] );
double *triangle_to_ref ( double tvert1[], double tvert2[], double tvert3[], 
  double t[] );
double *triangle_to_simplex ( double tvert1[], double tvert2[], double tvert3[], 
  double t[] );
void trianmap ( int numnodes, double vert1[], double vert2[], double vert3[], 
  double rnodes[], double whts[] );
void triasimp ( double x, double y, double &uout, double &vout );
void triasymq ( int n, double vert1[], double vert2[], double vert3[], 
  double rnodes[], double weights[], int numnodes );
void triasymq_gnuplot ( double vert1[], double vert2[], double vert3[], 
  int numnodes, double rnodes[], string header );
