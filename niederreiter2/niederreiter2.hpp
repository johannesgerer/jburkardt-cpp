//  NIEDERREITER2.H
//
//  Local Parameters:
//
//    int MAXDEG, the highest degree of polynomial
//    to be handled. 
//
//    int DIM_MAX, the maximum dimension that will be used.
//
//    int NBITS, the number of bits (not counting the sign) in a
//    fixed-point integer.
//
//    double RECIP, the multiplier which changes the
//    integers in NEXTQ into the required real values in QUASI.
//
# define MAXDEG 50
# define DIM_MAX 20
# define NBITS 31

//
//  Thanks to Bradley Keister for pointing out that the compiler might
//  have trouble using POW since it's not declared yet...
//
//static double RECIP = pow ( 2.0, -NBITS );

static double RECIP = 1.0 / ( double ) ( 1 << NBITS );

void calcc2 ( int dimen, int cj[DIM_MAX][NBITS] );

void calcv2 ( int maxv, int px_deg, int px[MAXDEG+1], int add[2][2], 
  int mul[2][2], int sub[2][2], int *b_deg, int b[MAXDEG+1], 
  int v[] );

void niederreiter2 ( int dim, int *seed, double quasi[] );

double *niederreiter2_generate ( int dim_num, int n, int *seed );

void plymul2 ( int add[2][2], int mul[2][2], int pa_deg, 
  int pa[MAXDEG+1], int pb_deg, int pb[MAXDEG+1], 
  int *pc_deg, int pc[MAXDEG+1] );

void setfld2 ( int add[2][2], int mul[2][2], int sub[2][2] );

void timestamp ( );
