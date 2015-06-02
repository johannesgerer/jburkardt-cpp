# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
int choose ( int n, int k );
void comp_next ( int n, int k, int a[], bool *more, int *h, int *t );
void hermite_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] );
void hermite_weights ( int order, double weight[] );
int i4_log_2 ( int i );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
string i4_to_string ( int i4, string format );
int i4vec_product ( int n, int a[] );
int *index_level_herm ( int level, int level_max, int dim_num, int point_num, 
  int grid_index[], int grid_base[] );
void level_to_order_open ( int dim_num, int level[], int order[] );
int *multigrid_index_z ( int dim_num, int order_1d[], int order_nd );
double *product_weight_herm ( int dim_num, int order_1d[], int order_nd );
double r8_epsilon ( );
double r8_huge ( );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] );
void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title );
int s_len_trim ( string s );
void sparse_grid_herm ( int dim_num, int level_max, int point_num, 
  double grid_weight[], double grid_point[] );
void sparse_grid_herm_index ( int dim_num, int level_max, int point_num, 
  int grid_index [], int grid_base[] );
int sparse_grid_herm_size ( int dim_num, int level_max );
void timestamp ( );
void vec_colex_next2 ( int dim_num, int base[], int a[], bool *more );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPARSE_GRID_HERMITE_DATASET.
//
//  Discussion:
//
//    This program computes a sparse grid quadrature rule based on a 1D
//    Gauss-Hermite rule and writes it to a file.. 
//
//    The user specifies:
//    * the spatial dimension of the quadrature region,
//    * the level that defines the Smolyak grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
{
  int dim;
  int dim_num;
  int level_max;
  int level_min;
  double pi = 3.141592653589793;
  int point;
  int point_num;
  double *r;
  string r_filename;
  double *w;
  string w_filename;
  double weight_sum;
  double *x;
  string x_filename;

  timestamp ( );
  cout << "\n";
  cout << "SPARSE_GRID_HERMITE_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Compute the abscissas and weights of a quadrature rule\n";
  cout << "  associated with a sparse grid derived from a Smolyak\n";
  cout << "  construction based on a 1D Gauss-Hermite rule.\n";
  cout << "\n";
  cout << "  Inputs to the program include:\n";
  cout << "\n";
  cout << "    DIM_NUM, the spatial dimension.\n";
  cout << "    (typically in the range of 2 to 10)\n";
  cout << "\n";
  cout << "    LEVEL_MAX, the level of the sparse grid.\n";
  cout << "    (typically in the range of 0, 1, 2, 3, ...\n";
  cout << "\n";
  cout << "  Output from the program includes:\n";
  cout << "\n";
  cout << "    * A printed table of the abscissas and weights.\n";
  cout << "\n";
  cout << "    * A set of 3 files that define the quadrature rule.\n";
  cout << "\n";
  cout << "    (1) herm_d?_level?_r.txt, the ranges;\n";
  cout << "    (2) herm_d?_level?_w.txt, the weights;\n";
  cout << "    (3) herm_d?_level?_x.txt, the abscissas.\n";
//
//  Get the spatial dimension.
//
  if ( 1 < argc )
  {
    dim_num = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of DIM_NUM (1 or greater)\n";
    cin >> dim_num;
  }

  cout << "\n";
  cout << "  Spatial dimension requested is = " << dim_num << "\n";
//
//  Get the level.
//
  if ( 2 < argc )
  {
    level_max = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of LEVEL_MAX (0 or greater).\n";
    cin >> level_max;
  }

  level_min = i4_max ( 0, level_max + 1 - dim_num );

  cout << "\n";
  cout << "  LEVEL_MIN is = " << level_min << "\n";
  cout << "  LEVEL_MAX is = " << level_max << "\n";
// 
//  How many distinct points will there be?
//
  point_num = sparse_grid_herm_size ( dim_num, level_max );

  cout << "\n";
  cout << "  The number of distinct abscissas in the\n";
  cout << "  quadrature rule is determined from the spatial\n";
  cout << "  dimension DIM_NUM and the level LEVEL_MAX.\n";
  cout << "  For the given input, this value will be = " << point_num << "\n";
//
//  Allocate memory.
//
  r = new double[dim_num*2];
  w = new double[point_num];
  x = new double[dim_num*point_num];
//
//  Compute the weights and points.
//
  for ( dim = 0; dim < dim_num; dim++ )
  {
    r[dim+0*dim_num] = - r8_huge ( );
    r[dim+1*dim_num] = + r8_huge ( );
  }

  sparse_grid_herm ( dim_num, level_max, point_num, w, x );

  r8mat_transpose_print_some ( dim_num, point_num, x, 1, 1, dim_num, 
    10, "  First 10 grid points:" );

  r8vec_print_some ( point_num, w, 1, 10, "  First 10 grid weights:" );

  weight_sum = 0.0;
  for ( point = 0; point < point_num; point++ )
  {
    weight_sum = weight_sum + w[point];
  }

  cout << "\n";
  cout << "  Weights sum to   " << weight_sum << "\n";
  cout << "  Correct value is " << sqrt ( pow ( pi, dim_num ) ) << "\n";
//
//  Construct appropriate file names.
//
  r_filename = "herm_d" + i4_to_string ( dim_num, "%d" )
    + "_level" + i4_to_string ( level_max, "%d" ) + "_r.txt";
  w_filename = "herm_d" + i4_to_string ( dim_num, "%d" )
    + "_level" + i4_to_string ( level_max, "%d" ) + "_w.txt";
  x_filename = "herm_d" + i4_to_string ( dim_num, "%d" )
    + "_level" + i4_to_string ( level_max, "%d" ) + "_x.txt";
//
//  Write the rule to files.
//
  cout << "\n";
  cout << "  Creating R file = \"" << r_filename << "\".\n";

  r8mat_write ( r_filename, dim_num, 2, r );

  cout << "  Creating W file = \"" << w_filename << "\".\n";

  r8mat_write ( w_filename, 1, point_num, w );
 
  cout << "  Creating X file = \"" << x_filename << "\".\n";

  r8mat_write ( x_filename, dim_num, point_num, x );

  delete [] r;
  delete [] w;
  delete [] x;

  cout << "\n";
  cout << "SPARSE_GRID_HERMITE_DATASET\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    CHOOSE computes the binomial coefficient C(N,K).
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in integer arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, are the values of N and K.
//
//    Output, int CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0;
  }
  else if ( mn == 0 )
  {
    value = 1;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( mx + i ) ) / i;
    }
  }

  return value;
}
//****************************************************************************80

void comp_next ( int n, int k, int a[], bool *more, int *h, int *t )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT computes the compositions of the integer N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is
//    3+2+1, another would be 6+0+0.
//
//    On the first call to this routine, set MORE = FALSE.  The routine
//    will compute the first element in the sequence of compositions, and
//    return it, as well as setting MORE = TRUE.  If more compositions
//    are desired, call again, and again.  Each time, the routine will
//    return with a new composition.
//
//    However, when the LAST composition in the sequence is computed 
//    and returned, the routine will reset MORE to FALSE, signaling that
//    the end of the sequence has been reached.
//
//    This routine originally used a SAVE statement to maintain the
//    variables H and T.  I have decided that it is safer
//    to pass these variables as arguments, even though the user should
//    never alter them.  This allows this routine to safely shuffle
//    between several ongoing calculations.
//
//
//    There are 28 compositions of 6 into three parts.  This routine will
//    produce those compositions in the following order:
//
//     I         A
//     -     ---------
//     1     6   0   0
//     2     5   1   0
//     3     4   2   0
//     4     3   3   0
//     5     2   4   0
//     6     1   5   0
//     7     0   6   0
//     8     5   0   1
//     9     4   1   1
//    10     3   2   1
//    11     2   3   1
//    12     1   4   1
//    13     0   5   1
//    14     4   0   2
//    15     3   1   2
//    16     2   2   2
//    17     1   3   2
//    18     0   4   2
//    19     3   0   3
//    20     2   1   3
//    21     1   2   3
//    22     0   3   3
//    23     2   0   4
//    24     1   1   4
//    25     0   2   4
//    26     1   0   5
//    27     0   1   5
//    28     0   0   6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2008
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool *MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
//    Input/output, int *H, *T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;

  if ( !( *more ) )
  {
    *t = n;
    *h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < *t )
    {
      *h = 0;
    }
    *h = *h + 1;
    *t = a[*h-1];
    a[*h-1] = 0;
    a[0] = *t - 1;
    a[*h] = a[*h] + 1;
  }

  *more = ( a[k-1] != n );

  return;
}
//****************************************************************************80

void hermite_abscissa ( int dim_num, int point_num, int grid_index[], 
  int grid_base[], double grid_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_ABSCISSA sets abscissas for multidimensional Gauss-Hermite quadrature.
//
//  Discussion:
//
//    The "nesting" as it occurs for Gauss-Hermite sparse grids simply
//    involves the use of a specified set of permissible orders for the
//    rule.  
//
//    The X array lists the (complete) Gauss-Legendre abscissas for rules 
//    of order 1, 3, 7, 15, 31, 63 or 127, in order. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], the index of the abscissa
//    from the rule, for each dimension and point.
//
//    Input, int GRID_BASE[DIM_NUM], the number of points used in the 
//    rule for a given dimension.
//
//    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the grid points of abscissas.
//
{
  int dim;
  int level;
  int point;
  int pointer;
  int skip[8] = { 0, 1, 4, 11, 26, 57, 120, 247 };
  double x[247] = {
    0.0E+00, 
   -0.122474487139158904909864203735E+01, 
    0.0E+00, 
    0.122474487139158904909864203735E+01, 
   -0.265196135683523349244708200652E+01, 
   -0.167355162876747144503180139830E+01, 
   -0.816287882858964663038710959027E+00, 
    0.0E+00, 
    0.816287882858964663038710959027E+00, 
    0.167355162876747144503180139830E+01, 
    0.265196135683523349244708200652E+01, 
   -0.449999070730939155366438053053E+01, 
   -0.366995037340445253472922383312E+01, 
   -0.296716692790560324848896036355E+01, 
   -0.232573248617385774545404479449E+01, 
   -0.171999257518648893241583152515E+01, 
   -0.113611558521092066631913490556E+01, 
   -0.565069583255575748526020337198E+00, 
    0.0E+00, 
    0.565069583255575748526020337198E+00, 
    0.113611558521092066631913490556E+01, 
    0.171999257518648893241583152515E+01, 
    0.232573248617385774545404479449E+01, 
    0.296716692790560324848896036355E+01, 
    0.366995037340445253472922383312E+01, 
    0.449999070730939155366438053053E+01, 
   -6.9956801237185402753248521473232E+00, 
   -6.2750787049428601427036567812530E+00, 
   -5.6739614446185883296332558789276E+00, 
   -5.1335955771123807045862968913996E+00, 
   -4.6315595063128599420667997654336E+00, 
   -4.1562717558181451724831352315314E+00, 
   -3.7007434032314694224497164589673E+00, 
   -3.2603207323135408104645401509648E+00, 
   -2.8316804533902054557015640151425E+00, 
   -2.4123177054804201051740184582119E+00, 
   -2.0002585489356389657975562598571E+00, 
   -1.5938858604721398261388419455550E+00, 
   -1.1918269983500464260821358649242E+00, 
   -0.79287697691530893968593032998830E+00, 
   -0.39594273647142311094670041663436E+00, 
    0.0000000000000000000000000000000E+00, 
    0.39594273647142311094670041663436E+00, 
    0.79287697691530893968593032998830E+00, 
    1.1918269983500464260821358649242E+00, 
    1.5938858604721398261388419455550E+00, 
    2.0002585489356389657975562598571E+00, 
    2.4123177054804201051740184582119E+00, 
    2.8316804533902054557015640151425E+00, 
    3.2603207323135408104645401509648E+00, 
    3.7007434032314694224497164589673E+00, 
    4.1562717558181451724831352315314E+00, 
    4.6315595063128599420667997654336E+00, 
    5.1335955771123807045862968913996E+00, 
    5.6739614446185883296332558789276E+00, 
    6.2750787049428601427036567812530E+00, 
    6.9956801237185402753248521473232E+00, 
   -10.435499877854168053468115427285E+00, 
   -9.8028759912974963635223935286507E+00, 
   -9.2792019543050391319404745506496E+00, 
   -8.8118581437284546442526628275570E+00, 
   -8.3807683451863219343010651043788E+00, 
   -7.9755950801420373181541806298501E+00, 
   -7.5901395198641066762479783194468E+00, 
   -7.2203167078889678461161324222529E+00, 
   -6.8632544331795368527353285876066E+00, 
   -6.5168348106821160605273395854042E+00, 
   -6.1794379922705969862418461787263E+00, 
   -5.8497884000810673462526582961482E+00, 
   -5.5268572526403031425047575122840E+00, 
   -5.2097979830408354861575136416263E+00, 
   -4.8979018644975742350745099214868E+00, 
   -4.5905665744435190229271294569091E+00, 
   -4.2872733352824404031727616199454E+00, 
   -3.9875699104197157485227052068068E+00, 
   -3.6910577000963465117322810559754E+00, 
   -3.3973817713303911852755941806287E+00, 
   -3.1062230279282566329138616746036E+00, 
   -2.8172919672837977750747135657355E+00, 
   -2.5303236304712010926855221718499E+00, 
   -2.2450734604812066298995918179330E+00, 
   -1.9613138583081485293922008411321E+00, 
   -1.6788312791720137520802800622638E+00, 
   -1.3974237486049625107570752063702E+00, 
   -1.1168987050996462690510970277840E+00, 
   -0.83707109558947615977737795461293E+00, 
   -0.55776166427908221668763665253822E+00, 
   -0.27879538567115223986687628627202E+00, 
    0.00000000000000000000000000000000E+00, 
    0.27879538567115223986687628627202E+00, 
    0.55776166427908221668763665253822E+00, 
    0.83707109558947615977737795461293E+00, 
    1.1168987050996462690510970277840E+00, 
    1.3974237486049625107570752063702E+00, 
    1.6788312791720137520802800622638E+00, 
    1.9613138583081485293922008411321E+00, 
    2.2450734604812066298995918179330E+00, 
    2.5303236304712010926855221718499E+00, 
    2.8172919672837977750747135657355E+00, 
    3.1062230279282566329138616746036E+00, 
    3.3973817713303911852755941806287E+00, 
    3.6910577000963465117322810559754E+00, 
    3.9875699104197157485227052068068E+00, 
    4.2872733352824404031727616199454E+00, 
    4.5905665744435190229271294569091E+00, 
    4.8979018644975742350745099214868E+00, 
    5.2097979830408354861575136416263E+00, 
    5.5268572526403031425047575122840E+00, 
    5.8497884000810673462526582961482E+00, 
    6.1794379922705969862418461787263E+00, 
    6.5168348106821160605273395854042E+00, 
    6.8632544331795368527353285876066E+00, 
    7.2203167078889678461161324222529E+00, 
    7.5901395198641066762479783194468E+00, 
    7.9755950801420373181541806298501E+00, 
    8.3807683451863219343010651043788E+00, 
    8.8118581437284546442526628275570E+00, 
    9.2792019543050391319404745506496E+00, 
    9.8028759912974963635223935286507E+00, 
    10.435499877854168053468115427285E+00, 
   -15.228338148167350978246954433464E+00, 
   -14.669595158833972632746354112896E+00, 
   -14.209085995284870755168244250887E+00, 
   -13.799722290211676634645246746673E+00, 
   -13.423518590070950062438258321855E+00, 
   -13.071208660474601901583995439649E+00, 
   -12.737235652415686338138003924072E+00, 
   -12.417939378869715805445879624069E+00, 
   -12.110749020947747600132123508132E+00, 
   -11.813772198267727195134584136191E+00, 
   -11.525565112572696599167888588564E+00, 
   -11.244994583785543445194384194300E+00, 
   -10.971150569840247423423040263881E+00, 
   -10.703288201027481347670940744690E+00, 
   -10.440787957772772867742591798027E+00, 
   -10.183127473450343888624126450357E+00, 
   -9.9298610495114250736847004273684E+00, 
   -9.6806044412474728038150712732737E+00, 
   -9.4350233389881650135019598506287E+00, 
   -9.1928244988460305715774195052527E+00, 
   -8.9537488108565404323807890169970E+00, 
   -8.7175658087076307363833999548548E+00, 
   -8.4840692689832473326097180339984E+00, 
   -8.2530736454457156579694124243888E+00, 
   -8.0244111514703375578594739796798E+00, 
   -7.7979293513870105420829120455591E+00, 
   -7.5734891556083454022834960763301E+00, 
   -7.3509631392269052701961258043733E+00, 
   -7.1302341220350710668064025713431E+00, 
   -6.9111939615465713197465633109366E+00, 
   -6.6937425208758294190074417381666E+00, 
   -6.4777867811645365448144903821487E+00, 
   -6.2632400742737354345609723857092E+00, 
   -6.0500214161419845694465474482388E+00, 
   -5.8380549248774187386601690807757E+00, 
   -5.6272693105464816659423455794909E+00, 
   -5.4175974259243240722848425872924E+00, 
   -5.2089758693153983587570258372239E+00, 
   -5.0013446320386360038520809107373E+00, 
   -4.7946467843764925009748509930857E+00, 
   -4.5888281947698372951606485031212E+00, 
   -4.3838372778464736294253744407459E+00, 
   -4.1796247675352031349421189892408E+00, 
   -3.9761435120673355916035814195920E+00, 
   -3.7733482881250526721004678400057E+00, 
   -3.5711956317782180447199756485249E+00, 
   -3.3696436841717397896643629240035E+00, 
   -3.1686520501953630191857798261495E+00, 
   -2.9681816685955910267761649521505E+00, 
   -2.7681946921824058801226545958892E+00, 
   -2.5686543769473501723144013022363E+00, 
   -2.3695249790490401080012474645702E+00, 
   -2.1707716587411506879498498083695E+00, 
   -1.9723603904195020079324743227565E+00, 
   -1.7742578780516791584676442103681E+00, 
   -1.5764314753267801315519597621879E+00, 
   -1.3788491099261778091441557053728E+00, 
   -1.1814792113700685848678583598423E+00, 
   -0.98429064194027277726568984213773E+00, 
   -0.78725263021825034151596831878971E+00, 
   -0.59033470680942102142230439346102E+00, 
   -0.39350664185130136568037826200185E+00, 
   -0.19673838392423251964272239737078E+00, 
    0.0000000000000000000000000000000E+00, 
    0.19673838392423251964272239737078E+00, 
    0.39350664185130136568037826200185E+00, 
    0.59033470680942102142230439346102E+00, 
    0.78725263021825034151596831878971E+00, 
    0.98429064194027277726568984213773E+00, 
    1.1814792113700685848678583598423E+00, 
    1.3788491099261778091441557053728E+00, 
    1.5764314753267801315519597621879E+00, 
    1.7742578780516791584676442103681E+00, 
    1.9723603904195020079324743227565E+00, 
    2.1707716587411506879498498083695E+00, 
    2.3695249790490401080012474645702E+00, 
    2.5686543769473501723144013022363E+00, 
    2.7681946921824058801226545958892E+00, 
    2.9681816685955910267761649521505E+00, 
    3.1686520501953630191857798261495E+00, 
    3.3696436841717397896643629240035E+00, 
    3.5711956317782180447199756485249E+00, 
    3.7733482881250526721004678400057E+00, 
    3.9761435120673355916035814195920E+00, 
    4.1796247675352031349421189892408E+00, 
    4.3838372778464736294253744407459E+00, 
    4.5888281947698372951606485031212E+00, 
    4.7946467843764925009748509930857E+00, 
    5.0013446320386360038520809107373E+00, 
    5.2089758693153983587570258372239E+00, 
    5.4175974259243240722848425872924E+00, 
    5.6272693105464816659423455794909E+00, 
    5.8380549248774187386601690807757E+00, 
    6.0500214161419845694465474482388E+00, 
    6.2632400742737354345609723857092E+00, 
    6.4777867811645365448144903821487E+00, 
    6.6937425208758294190074417381666E+00, 
    6.9111939615465713197465633109366E+00, 
    7.1302341220350710668064025713431E+00, 
    7.3509631392269052701961258043733E+00, 
    7.5734891556083454022834960763301E+00, 
    7.7979293513870105420829120455591E+00, 
    8.0244111514703375578594739796798E+00, 
    8.2530736454457156579694124243888E+00, 
    8.4840692689832473326097180339984E+00, 
    8.7175658087076307363833999548548E+00, 
    8.9537488108565404323807890169970E+00, 
    9.1928244988460305715774195052527E+00, 
    9.4350233389881650135019598506287E+00, 
    9.6806044412474728038150712732737E+00, 
    9.9298610495114250736847004273684E+00, 
    10.183127473450343888624126450357E+00, 
    10.440787957772772867742591798027E+00, 
    10.703288201027481347670940744690E+00, 
    10.971150569840247423423040263881E+00, 
    11.244994583785543445194384194300E+00, 
    11.525565112572696599167888588564E+00, 
    11.813772198267727195134584136191E+00, 
    12.110749020947747600132123508132E+00, 
    12.417939378869715805445879624069E+00, 
    12.737235652415686338138003924072E+00, 
    13.071208660474601901583995439649E+00, 
    13.423518590070950062438258321855E+00, 
    13.799722290211676634645246746673E+00, 
    14.209085995284870755168244250887E+00, 
    14.669595158833972632746354112896E+00, 
    15.228338148167350978246954433464E+00  
  };

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( grid_base[dim] < 0 )
    {
      cout << "\n";
      cout << "HERMITE_ABSCISSA - Fatal error!\n";
      cout << "  Some base values are less than 0.\n";
      exit ( 1 );
    }
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( 63 < grid_base[dim] )
    {
      cout << "\n";
      cout << "HERMITE_ABSCISSA - Fatal error!\n";
      cout << "  Some base values are greater than 63.\n";
      exit ( 1 );
    }
  }

  for ( point = 0; point < point_num; point++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      level = i4_log_2 ( grid_base[dim] + 1 );

      pointer = skip[level] + ( grid_index[dim+point*dim_num] + grid_base[dim] );

      grid_point[dim+point*dim_num] = x[pointer];
    }
  }

  return;
}
//****************************************************************************80

void hermite_weights ( int order, double weight[] )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_WEIGHTS returns weights for certain Gauss-Hermite quadrature rules.
//
//  Discussion:
//
//    The allowed orders are 1, 3, 7, 15, 31, 63 and 127.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
//
//    Output, double WEIGHT[ORDER], the weights.
//    The weights are positive, symmetric and should sum to SQRT(PI).
//
{
  if ( order == 1 )
  {
    weight[1-1] = 1.77245385090551602729816748334E+00;
  }
  else if ( order == 3 )
  {
    weight[1-1] = 0.295408975150919337883027913890E+00;
    weight[2-1] = 0.118163590060367735153211165556E+01;
    weight[3-1] = 0.295408975150919337883027913890E+00;
  }
  else if ( order == 7 )
  {
    weight[1-1] = 0.971781245099519154149424255939E-03;
    weight[2-1] = 0.545155828191270305921785688417E-01;
    weight[3-1] = 0.425607252610127800520317466666E+00;
    weight[4-1] = 0.810264617556807326764876563813E+00;
    weight[5-1] = 0.425607252610127800520317466666E+00;
    weight[6-1] = 0.545155828191270305921785688417E-01;
    weight[7-1] = 0.971781245099519154149424255939E-03;
  }
  else if ( order == 15 )
  {
    weight[1-1] =  0.152247580425351702016062666965E-08;
    weight[2-1] =  0.105911554771106663577520791055E-05;
    weight[3-1] =  0.100004441232499868127296736177E-03;
    weight[4-1] =  0.277806884291277589607887049229E-02;
    weight[5-1] =  0.307800338725460822286814158758E-01;
    weight[6-1] =  0.158488915795935746883839384960E+00;
    weight[7-1] =  0.412028687498898627025891079568E+00;
    weight[8-1] =  0.564100308726417532852625797340E+00;
    weight[9-1] =  0.412028687498898627025891079568E+00;
    weight[10-1] = 0.158488915795935746883839384960E+00;
    weight[11-1] = 0.307800338725460822286814158758E-01;
    weight[12-1] = 0.277806884291277589607887049229E-02;
    weight[13-1] = 0.100004441232499868127296736177E-03;
    weight[14-1] = 0.105911554771106663577520791055E-05;
    weight[15-1] = 0.152247580425351702016062666965E-08;
  }
  else if ( order == 31 )
  {
    weight[  1-1] =   0.46189683944498305857470556847735E-21;
    weight[  2-1] =   0.51106090079112519643027197715274E-17;
    weight[  3-1] =   0.58995564987355133075257722133966E-14;
    weight[  4-1] =   0.18603735214463569590294465062239E-11;
    weight[  5-1] =   0.23524920032013205739850619940094E-09;
    weight[  6-1] =   0.14611988344865057576066495091513E-07;
    weight[  7-1] =   0.50437125589241034841778074689627E-06;
    weight[  8-1] =   0.10498602757642934202945441341697E-04;
    weight[  9-1] =   0.13952090395003623854995664958146E-03;
    weight[ 10-1] =   0.12336833073030489880608311394968E-02;
    weight[ 11-1] =   0.74827999140119116765002499116934E-02;
    weight[ 12-1] =   0.31847230731201222775249585776902E-01;
    weight[ 13-1] =   0.96717948160569462991143316029341E-01;
    weight[ 14-1] =   0.21213278866810461318136114862419E+00;
    weight[ 15-1] =   0.33877265789305344906000174083214E+00;
    weight[ 16-1] =   0.39577855609737786462923720809676E+00;
    weight[ 17-1] =   0.33877265789305344906000174083214E+00;
    weight[ 18-1] =   0.21213278866810461318136114862419E+00;
    weight[ 19-1] =   0.96717948160569462991143316029341E-01;
    weight[ 20-1] =   0.31847230731201222775249585776902E-01;
    weight[ 21-1] =   0.74827999140119116765002499116934E-02;
    weight[ 22-1] =   0.12336833073030489880608311394968E-02;
    weight[ 23-1] =   0.13952090395003623854995664958146E-03;
    weight[ 24-1] =   0.10498602757642934202945441341697E-04;
    weight[ 25-1] =   0.50437125589241034841778074689627E-06;
    weight[ 26-1] =   0.14611988344865057576066495091513E-07;
    weight[ 27-1] =   0.23524920032013205739850619940094E-09;
    weight[ 28-1] =   0.18603735214463569590294465062239E-11;
    weight[ 29-1] =   0.58995564987355133075257722133966E-14;
    weight[ 30-1] =   0.51106090079112519643027197715274E-17;
    weight[ 31-1] =   0.46189683944498305857470556847735E-21;
  }
  else if ( order == 63 )
  {
    weight[  1-1] =   0.37099206434787551197827130470031E-47;
    weight[  2-1] =   0.10400778615192299534481914814892E-41;
    weight[  3-1] =   0.19796804708258311251124226474396E-37;
    weight[  4-1] =   0.84687478191640015120141181138947E-34;
    weight[  5-1] =   0.13071305930779945903630127634063E-30;
    weight[  6-1] =   0.93437837175367456929765381518998E-28;
    weight[  7-1] =   0.36027426635173044862245783257252E-25;
    weight[  8-1] =   0.82963863115951789374753323156164E-23;
    weight[  9-1] =   0.12266629909105281472971700203949E-20;
    weight[ 10-1] =   0.12288435628797061539461585325494E-18;
    weight[ 11-1] =   0.86925536958188009075932426691516E-17;
    weight[ 12-1] =   0.44857058689176221240330804981619E-15;
    weight[ 13-1] =   0.17335817955735154599902643794700E-13;
    weight[ 14-1] =   0.51265062385038307838565047455223E-12;
    weight[ 15-1] =   0.11808921844532942490513037158404E-10;
    weight[ 16-1] =   0.21508698297808025739828859845140E-09;
    weight[ 17-1] =   0.31371929535285447801497640621672E-08;
    weight[ 18-1] =   0.37041625984781705796752840204084E-07;
    weight[ 19-1] =   0.35734732949879669663960738150956E-06;
    weight[ 20-1] =   0.28393114498380927832990899215541E-05;
    weight[ 21-1] =   0.18709113003730498008961134765721E-04;
    weight[ 22-1] =   0.10284880800653635546698378640623E-03;
    weight[ 23-1] =   0.47411702610173128107201781718693E-03;
    weight[ 24-1] =   0.18409222622384813438539657470055E-02;
    weight[ 25-1] =   0.60436044551187631655712178246467E-02;
    weight[ 26-1] =   0.16829299199599730926458559757600E-01;
    weight[ 27-1] =   0.39858264027692992170237391875317E-01;
    weight[ 28-1] =   0.80467087993950415219587554532823E-01;
    weight[ 29-1] =   0.13871950817615293377792092082674E+00;
    weight[ 30-1] =   0.20448695346833761570957197160475E+00;
    weight[ 31-1] =   0.25799889943058042204920467417642E+00;
    weight[ 32-1] =   0.27876694884838411919175686949858E+00;
    weight[ 33-1] =   0.25799889943058042204920467417642E+00;
    weight[ 34-1] =   0.20448695346833761570957197160475E+00;
    weight[ 35-1] =   0.13871950817615293377792092082674E+00;
    weight[ 36-1] =   0.80467087993950415219587554532823E-01;
    weight[ 37-1] =   0.39858264027692992170237391875317E-01;
    weight[ 38-1] =   0.16829299199599730926458559757600E-01;
    weight[ 39-1] =   0.60436044551187631655712178246467E-02;
    weight[ 40-1] =   0.18409222622384813438539657470055E-02;
    weight[ 41-1] =   0.47411702610173128107201781718693E-03;
    weight[ 42-1] =   0.10284880800653635546698378640623E-03;
    weight[ 43-1] =   0.18709113003730498008961134765721E-04;
    weight[ 44-1] =   0.28393114498380927832990899215541E-05;
    weight[ 45-1] =   0.35734732949879669663960738150956E-06;
    weight[ 46-1] =   0.37041625984781705796752840204084E-07;
    weight[ 47-1] =   0.31371929535285447801497640621672E-08;
    weight[ 48-1] =   0.21508698297808025739828859845140E-09;
    weight[ 49-1] =   0.11808921844532942490513037158404E-10;
    weight[ 50-1] =   0.51265062385038307838565047455223E-12;
    weight[ 51-1] =   0.17335817955735154599902643794700E-13;
    weight[ 52-1] =   0.44857058689176221240330804981619E-15;
    weight[ 53-1] =   0.86925536958188009075932426691516E-17;
    weight[ 54-1] =   0.12288435628797061539461585325494E-18;
    weight[ 55-1] =   0.12266629909105281472971700203949E-20;
    weight[ 56-1] =   0.82963863115951789374753323156164E-23;
    weight[ 57-1] =   0.36027426635173044862245783257252E-25;
    weight[ 58-1] =   0.93437837175367456929765381518998E-28;
    weight[ 59-1] =   0.13071305930779945903630127634063E-30;
    weight[ 60-1] =   0.84687478191640015120141181138947E-34;
    weight[ 61-1] =   0.19796804708258311251124226474396E-37;
    weight[ 62-1] =   0.10400778615192299534481914814892E-41;
    weight[ 63-1] =   0.37099206434787551197827130470031E-47;
  }
  else if ( order == 127 )
  {
    weight[  1-1] =   0.12504497577050595552677230002883E-100;
    weight[  2-1] =   0.17272798059419131415318615789672E-93;
    weight[  3-1] =   0.89321681571986548608031150791499E-88;
    weight[  4-1] =   0.77306185240893578449625186483810E-83;
    weight[  5-1] =   0.20143957652648255497735460506196E-78;
    weight[  6-1] =   0.21503714733610239701351039429345E-74;
    weight[  7-1] =   0.11341924208594594813715533569504E-70;
    weight[  8-1] =   0.33489139011795051950683388483136E-67;
    weight[  9-1] =   0.60486548964016681064424451668405E-64;
    weight[ 10-1] =   0.71375092946352177824971347343892E-61;
    weight[ 11-1] =   0.57884563374885556636801095624030E-58;
    weight[ 12-1] =   0.33581166223858230300409326551248E-55;
    weight[ 13-1] =   0.14394641949253923568603163698953E-52;
    weight[ 14-1] =   0.46821808383216117724080263903889E-50;
    weight[ 15-1] =   0.11817054440684264071348471955361E-47;
    weight[ 16-1] =   0.23581659156008927203181682045005E-45;
    weight[ 17-1] =   0.37814427940797540210712758405540E-43;
    weight[ 18-1] =   0.49411031115771638145610738414006E-41;
    weight[ 19-1] =   0.53255303775425059266087298458297E-39;
    weight[ 20-1] =   0.47854390680131484999315199332765E-37;
    weight[ 21-1] =   0.36191883445952356128627543209554E-35;
    weight[ 22-1] =   0.23232083386343554805352497446119E-33;
    weight[ 23-1] =   0.12753331411008716683688974281454E-31;
    weight[ 24-1] =   0.60277753850758742112436095241270E-30;
    weight[ 25-1] =   0.24679773241777200207460855084439E-28;
    weight[ 26-1] =   0.88019567691698482573264198727415E-27;
    weight[ 27-1] =   0.27482489212040561315005725890593E-25;
    weight[ 28-1] =   0.75468218903085486125222816438456E-24;
    weight[ 29-1] =   0.18303134636280466270545996891835E-22;
    weight[ 30-1] =   0.39355990860860813085582448449811E-21;
    weight[ 31-1] =   0.75293161638581191068419292570042E-20;
    weight[ 32-1] =   0.12857997786722855037584105682618E-18;
    weight[ 33-1] =   0.19659326888445857792541925311450E-17;
    weight[ 34-1] =   0.26986511907214101894995783364250E-16;
    weight[ 35-1] =   0.33344414303198856330118301113874E-15;
    weight[ 36-1] =   0.37173303125150639885726463109574E-14;
    weight[ 37-1] =   0.37473954472839737091885387788983E-13;
    weight[ 38-1] =   0.34230094493397259538669512076007E-12;
    weight[ 39-1] =   0.28385303724993373166810860630552E-11;
    weight[ 40-1] =   0.21406920290454669208938772802828E-10;
    weight[ 41-1] =   0.14706331273431716244229273183839E-09;
    weight[ 42-1] =   0.92173940967434659264335883218167E-09;
    weight[ 43-1] =   0.52781663936972714041837056042506E-08;
    weight[ 44-1] =   0.27650497044951117835905283127679E-07;
    weight[ 45-1] =   0.13267855842539464770913063113371E-06;
    weight[ 46-1] =   0.58380944276113062188573331195042E-06;
    weight[ 47-1] =   0.23581561724775629112332165335800E-05;
    weight[ 48-1] =   0.87524468034280444703919485644809E-05;
    weight[ 49-1] =   0.29876790535909012274846532159647E-04;
    weight[ 50-1] =   0.93874435720072545206729594267039E-04;
    weight[ 51-1] =   0.27170762627931172053444716883938E-03;
    weight[ 52-1] =   0.72493929742498358979684249380921E-03;
    weight[ 53-1] =   0.17841208326763432884316727108264E-02;
    weight[ 54-1] =   0.40524855186046131499765636276283E-02;
    weight[ 55-1] =   0.85000263041544110385806705526917E-02;
    weight[ 56-1] =   0.16471142241609687824005585301760E-01;
    weight[ 57-1] =   0.29499296248213632269675010319119E-01;
    weight[ 58-1] =   0.48847387114300011006959603975676E-01;
    weight[ 59-1] =   0.74807989768583731416517226905270E-01;
    weight[ 60-1] =   0.10598520508090929403834368934301E+00;
    weight[ 61-1] =   0.13893945309051540832066283010510E+00;
    weight[ 62-1] =   0.16856236074207929740526975049765E+00;
    weight[ 63-1] =   0.18927849580120432177170145550076E+00;
    weight[ 64-1] =   0.19673340688823289786163676995151E+00;
    weight[ 65-1] =   0.18927849580120432177170145550076E+00;
    weight[ 66-1] =   0.16856236074207929740526975049765E+00;
    weight[ 67-1] =   0.13893945309051540832066283010510E+00;
    weight[ 68-1] =   0.10598520508090929403834368934301E+00;
    weight[ 69-1] =   0.74807989768583731416517226905270E-01;
    weight[ 70-1] =   0.48847387114300011006959603975676E-01;
    weight[ 71-1] =   0.29499296248213632269675010319119E-01;
    weight[ 72-1] =   0.16471142241609687824005585301760E-01;
    weight[ 73-1] =   0.85000263041544110385806705526917E-02;
    weight[ 74-1] =   0.40524855186046131499765636276283E-02;
    weight[ 75-1] =   0.17841208326763432884316727108264E-02;
    weight[ 76-1] =   0.72493929742498358979684249380921E-03;
    weight[ 77-1] =   0.27170762627931172053444716883938E-03;
    weight[ 78-1] =   0.93874435720072545206729594267039E-04;
    weight[ 79-1] =   0.29876790535909012274846532159647E-04;
    weight[ 80-1] =   0.87524468034280444703919485644809E-05;
    weight[ 81-1] =   0.23581561724775629112332165335800E-05;
    weight[ 82-1] =   0.58380944276113062188573331195042E-06;
    weight[ 83-1] =   0.13267855842539464770913063113371E-06;
    weight[ 84-1] =   0.27650497044951117835905283127679E-07;
    weight[ 85-1] =   0.52781663936972714041837056042506E-08;
    weight[ 86-1] =   0.92173940967434659264335883218167E-09;
    weight[ 87-1] =   0.14706331273431716244229273183839E-09;
    weight[ 88-1] =   0.21406920290454669208938772802828E-10;
    weight[ 89-1] =   0.28385303724993373166810860630552E-11;
    weight[ 90-1] =   0.34230094493397259538669512076007E-12;
    weight[ 91-1] =   0.37473954472839737091885387788983E-13;
    weight[ 92-1] =   0.37173303125150639885726463109574E-14;
    weight[ 93-1] =   0.33344414303198856330118301113874E-15;
    weight[ 94-1] =   0.26986511907214101894995783364250E-16;
    weight[ 95-1] =   0.19659326888445857792541925311450E-17;
    weight[ 96-1] =   0.12857997786722855037584105682618E-18;
    weight[ 97-1] =   0.75293161638581191068419292570042E-20;
    weight[ 98-1] =   0.39355990860860813085582448449811E-21;
    weight[ 99-1] =   0.18303134636280466270545996891835E-22;
    weight[100-1] =   0.75468218903085486125222816438456E-24;
    weight[101-1] =   0.27482489212040561315005725890593E-25;
    weight[102-1] =   0.88019567691698482573264198727415E-27;
    weight[103-1] =   0.24679773241777200207460855084439E-28;
    weight[104-1] =   0.60277753850758742112436095241270E-30;
    weight[105-1] =   0.12753331411008716683688974281454E-31;
    weight[106-1] =   0.23232083386343554805352497446119E-33;
    weight[107-1] =   0.36191883445952356128627543209554E-35;
    weight[108-1] =   0.47854390680131484999315199332765E-37;
    weight[109-1] =   0.53255303775425059266087298458297E-39;
    weight[110-1] =   0.49411031115771638145610738414006E-41;
    weight[111-1] =   0.37814427940797540210712758405540E-43;
    weight[112-1] =   0.23581659156008927203181682045005E-45;
    weight[113-1] =   0.11817054440684264071348471955361E-47;
    weight[114-1] =   0.46821808383216117724080263903889E-50;
    weight[115-1] =   0.14394641949253923568603163698953E-52;
    weight[116-1] =   0.33581166223858230300409326551248E-55;
    weight[117-1] =   0.57884563374885556636801095624030E-58;
    weight[118-1] =   0.71375092946352177824971347343892E-61;
    weight[119-1] =   0.60486548964016681064424451668405E-64;
    weight[120-1] =   0.33489139011795051950683388483136E-67;
    weight[121-1] =   0.11341924208594594813715533569504E-70;
    weight[122-1] =   0.21503714733610239701351039429345E-74;
    weight[123-1] =   0.20143957652648255497735460506196E-78;
    weight[124-1] =   0.77306185240893578449625186483810E-83;
    weight[125-1] =   0.89321681571986548608031150791499E-88;
    weight[126-1] =   0.17272798059419131415318615789672E-93;
    weight[127-1] =   0.12504497577050595552677230002883E-100;
  }
  else
  {
    cout << "\n";
    cout << "HERMITE_WEIGHTS - Fatal error!\n";
    cout << "  Illegal value of ORDER = " << order << "\n";
    cout << "  Legal values are 1, 3, 7, 15, 31, 63 and 127.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

int i4_log_2 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
//
//  Example:
//
//        I  I4_LOG_10
//    -----  --------
//        0    0
//        1    0
//        2    1
//        3    1
//        4    2
//        5    2
//        7    2
//        8    3
//        9    3
//     1000    9
//     1024   10
//
//  Discussion:
//
//    I4_LOG_2 ( I ) + 1 is the number of binary digits in I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number whose logarithm base 2 is desired.
//
//    Output, int I4_LOG_2, the integer part of the logarithm base 2 of
//    the absolute value of X.
//
{
  int i_abs;
  int two_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    two_pow = 2;

    i_abs = abs ( i );

    while ( two_pow <= i_abs )
    {
      value = value + 1;
      two_pow = two_pow * 2;
    }

  }

  return value;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MAX, the larger of i1 and i2.
//
{
  int value;

  if ( i2 < i1 ) 
  {
    value = i1;
  }
  else 
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
//
{
  int value;

  if ( i1 < i2 ) 
  {
    value = i1;
  }
  else 
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Formula:
//
//    If 
//      NREM = I4_MODP ( I, J ) 
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//  Comments:
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Examples:
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
// 
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is 
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cout << "\n";
    cout << "I4_MODP - Fatal error!\n";
    cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

string i4_to_string ( int i4, string format )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  char i4_char[80];
  string i4_string;

  sprintf ( i4_char, format.c_str ( ), i4 );

  i4_string = string ( i4_char );

  return i4_string;
}
//****************************************************************************80

int i4vec_product ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRODUCT multiplies the entries of an I4VEC.
//
//  Example:
//
//    A = ( 1, 2, 3, 4 )
//
//    I4VEC_PRODUCT = 24
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector
//
//    Output, int I4VEC_PRODUCT, the product of the entries of A.
//
{
  int i;
  int product;

  product = 1;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}
//****************************************************************************80

int *index_level_herm ( int level, int level_max, int dim_num, int point_num, 
  int grid_index[], int grid_base[] )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_LEVEL_HERM: determine first level at which given index is generated.
//
//  Discussion:
//
//    We are constructing a sparse grid of Gauss-Hermite points.  The grid
//    is built up of product grids, with a characteristic LEVEL.  
//
//    We are concerned with identifying points in this product grid which
//    have actually been generated previously, on a lower value of LEVEL.
//
//    This routine determines the lowest value of LEVEL at which each of
//    the input points would be generated.
//
//    In 1D, given LEVEL, the number of points is ORDER = 2**(LEVEL+1) + 1,
//    (except that LEVEL = 0 implies ORDER = 1//), the BASE is (ORDER-1)/2, 
//    and the point INDEX values range from -BASE to +BASE.
//
//    The values of INDEX and BASE allow us to determine the abstract
//    properties of the point.  In particular, if INDEX is 0, the corresponding
//    Gauss-Legendre abscissa is 0, the special "nested" value we need
//    to take care of.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int LEVEL, the level at which these points were 
//    generated.  LEVEL_MIN <= LEVEL <= LEVEL_MAX.
//
//    Input, int LEVEL_MAX, the maximum level.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points to be tested.
//
//    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], the indices of the 
//    points to be tested.
//
//    Input, int GRID_BASE[DIM_NUM], the "base", which is essentially
//    the denominator of the index.
//
//    Output, int INDEX_LEVEL_HERM[POINT_NUM], the value of LEVEL at 
//    which the point would first be generated.  This will be the same as
//    the input value of LEVEL, unless the point has an INDEX of 0 and
//    a corresponding BASE that is NOT zero.
//
{
  int dim;
  int *grid_level;
  int level_min;
  int point;
  
  grid_level = new int[point_num];
  
  level_min = i4_max ( 0, level_max + 1 - dim_num );
//
//  If a point has a DIM-th component whose INDEX is 0, then the 
//  value of LEVEL at which this point would first be generated is
//  less than LEVEL, unless the DIM-th component of GRID_BASE is 0.
//
  for ( point = 0; point < point_num; point++ )
  {
    grid_level[point] = i4_max ( level, level_min );

    for ( dim = 0; dim < dim_num; dim++ )
    {
      if ( grid_index[dim+point*dim_num] == 0 )
      {
        grid_level[point] = i4_max ( grid_level[point] - grid_base[dim], level_min );
      }
    }
  }

  return grid_level;
}
//****************************************************************************80

void level_to_order_open ( int dim_num, int level[], int order[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEVEL_TO_ORDER_OPEN converts a level to an order for open rules.
//
//  Discussion:
//
//    Sparse grids can naturally be nested.  A natural scheme is to use
//    a series of one-dimensional rules arranged in a series of "levels"
//    whose order roughly doubles with each step.
//
//    The arrangement described here works naturally for the Fejer Type 1,
//    Fejer Type 2, Newton Cotes Open, Newton Cotes Half Open,
//    and Gauss-Patterson rules.  It also can be used, partially, to describe
//    the growth of Gauss-Legendre and Gauss-Hermite rules.
//
//    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single 
//    point at the center, and for all values afterwards, we use the 
//    relationship
//
//      ORDER = 2**(LEVEL+1) - 1.
//
//    The following table shows how the growth will occur:
//
//    Level    Order
//
//    0          1
//    1          3 =  4 - 1
//    2          7 =  8 - 1
//    3         15 = 16 - 1
//    4         31 = 32 - 1
//    5         63 = 64 - 1
//
//    For the Fejer Type 1, Fejer Type 2, Newton Cotes Open, 
//    Newton Cotes Open Half, and Gauss-Patterson rules, the point growth is
//    nested.  If we have ORDER points on a particular LEVEL, the next level 
//    includes all these old points, plus ORDER+1 new points, formed in the 
//    gaps between successive pairs of old points plus an extra point at each 
//    end.
//
//    Level    Order = New + Old
//
//    0          1   =  1  +  0
//    1          3   =  2  +  1
//    2          7   =  4  +  3
//    3         15   =  8  +  7
//    4         31   = 16  + 15
//    5         63   = 32  + 31
//
//    If we use a series of Gauss-Legendre or Gauss-Hermite rules, then there 
//    is almost no nesting, except that the central point is shared.  If we 
//    insist on producing a comparable series of such points, then the 
//    "nesting" behavior is as follows:
//
//    Level    Order = New + Old
//
//    0          1   =  1  +  0
//    1          3   =  2  +  1
//    2          7   =  6  +  1
//    3         15   = 14  +  1
//    4         31   = 30  +  1
//    5         63   = 62  +  1
//
//    Moreover, if we consider ALL the points used in such a set of "nested" 
//    Gauss-Legendre or Gauss-Hermite rules, then we must sum the "NEW" column,
//    and we see that we get roughly twice as many points as for the truly 
//    nested rules.
//
//    In this routine, we assume that a vector of levels is given,
//    and the corresponding orders are desired.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL[DIM_NUM], the nesting level.
//
//    Output, int ORDER[DIM_NUM], the order (number of points) 
//    of the rule.
//
{
  int dim;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( level[dim] < 0 )
    {
      order[dim] = -1;
    }
    else if ( level[dim] == 0 )
    {
      order[dim] = 1;
    }
    else
    {
      order[dim] = i4_power ( 2, level[dim] + 1 ) - 1 ;
    }
  }
  return;
}
//****************************************************************************80

int *multigrid_index_z ( int dim_num, int order_1d[], int order_nd )

//****************************************************************************80
//
//  Purpose:
//
//    MULTIGRID_INDEX_Z returns an indexed multidimensional grid.
//
//  Discussion:
//
//    For dimension DIM, the number of points is ORDER_1D[DIM].
//
//    We assume that ORDER_1D[DIM] is an odd number,
//      ORDER_1D[DIM] = N = 2 * M + 1
//    so that the points have coordinates
//      -M/M, -(M-1)/M, ..., -1/M, 0/M, 1/M, 2/M, 3/M, ..., (M-1)/M, M/M.
//    and we index them as
//      -M,   -(M-1),        -1,   0,   1,   2,   3,   ...,  M-1,    M.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension of the points.
//
//    Input, int ORDER_1D[DIM_NUM], the order of the
//    rule in each dimension.
//
//    Input, int ORDER_ND, the product of the entries of ORDER_1D.
//
//    Output, int INDX[DIM_NUM*ORDER_ND], the indices of the points in
//    the grid.  The second dimension of this array is equal to the
//    product of the entries of ORDER_1D.
//
{
  int *a;
  int dim;
  bool more;
  int p;
  int *indx;

  indx = new int[dim_num*order_nd];
  a = new int[dim_num];
  more = false;
  p = 0;

  for ( ; ; )
  {
    vec_colex_next2 ( dim_num, order_1d, a, &more );

    if ( !more )
    {
      break;
    }
//
//  The values of A(DIM) are between 0 and ORDER_1D(DIM)-1 = N - 1 = 2 * M.
//  Subtracting M sets the range to -M to +M, as we wish.
//
    for ( dim = 0; dim < dim_num; dim++ )
    {
      indx[dim+p*dim_num] = a[dim] - ( order_1d[dim] - 1 ) / 2;
    }
    p = p + 1;
  }

  delete [] a;

  return indx;
}
//****************************************************************************80

double *product_weight_herm ( int dim_num, int order_1d[], int order_nd )

//****************************************************************************80
//
//  Purpose:
//
//    PRODUCT_WEIGHT_HERM: weights for a product Gauss-Hermite rule.
//
//  Discussion:
//
//    This routine computes the weights for a quadrature rule which is
//    a product of 1D Gauss-Hermite rules of varying order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int ORDER_1D[DIM_NUM], the order of the 1D rules.
//
//    Input, int ORDER_ND, the order of the product rule.
//
//    Output, double PRODUCT_WEIGHT_HERM[ORDER_ND], the product rule weights.
//
{
  int dim;
  int order;
  double *w_1d;
  double *w_nd;

  w_nd = new double[order_nd];
  
  for ( order = 0; order < order_nd; order++ )
  {
    w_nd[order] = 1.0;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    w_1d = new double[order_1d[dim]];
    
    hermite_weights ( order_1d[dim], w_1d );

    r8vec_direct_product2 ( dim, order_1d[dim], w_1d, dim_num, 
      order_nd, w_nd );

    delete [] w_1d;
  }
  return w_nd;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, an optional title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
//****************************************************************************80

void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    To explain what is going on here, suppose we had to construct
//    a multidimensional quadrature rule as the product of K rules
//    for 1D quadrature.
//
//    The product rule will be represented as a list of points and weights.
//
//    The J-th item in the product rule will be associated with
//      item J1 of 1D rule 1,
//      item J2 of 1D rule 2, 
//      ..., 
//      item JK of 1D rule K.
//
//    In particular, 
//      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
//    and
//      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
//
//    So we can construct the quadrature rule if we can properly
//    distribute the information in the 1D quadrature rules.
//
//    This routine carries out that task for the weights W.
//
//    Another way to do this would be to compute, one by one, the
//    set of all possible indices (J1,J2,...,JK), and then index
//    the appropriate information.  An advantage of the method shown
//    here is that you can process the K-th set of information and
//    then discard it.
//
//  Example:
//
//    Rule 1: 
//      Order = 4
//      W(1:4) = ( 2, 3, 5, 7 )
//
//    Rule 2:
//      Order = 3
//      W(1:3) = ( 11, 13, 17 )
//
//    Rule 3:
//      Order = 2
//      W(1:2) = ( 19, 23 )
//
//    Product Rule:
//      Order = 24
//      W(1:24) =
//        ( 2 * 11 * 19 )
//        ( 3 * 11 * 19 )
//        ( 4 * 11 * 19 )
//        ( 7 * 11 * 19 )
//        ( 2 * 13 * 19 )
//        ( 3 * 13 * 19 )
//        ( 5 * 13 * 19 )
//        ( 7 * 13 * 19 )
//        ( 2 * 17 * 19 )
//        ( 3 * 17 * 19 )
//        ( 5 * 17 * 19 )
//        ( 7 * 17 * 19 )
//        ( 2 * 11 * 23 )
//        ( 3 * 11 * 23 )
//        ( 5 * 11 * 23 )
//        ( 7 * 11 * 23 )
//        ( 2 * 13 * 23 )
//        ( 3 * 13 * 23 )
//        ( 5 * 13 * 23 )
//        ( 7 * 13 * 23 )
//        ( 2 * 17 * 23 )
//        ( 3 * 17 * 23 )
//        ( 5 * 17 * 23 )
//        ( 7 * 17 * 23 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int FACTOR_INDEX, the index of the factor being processed.
//    The first factor processed must be factor 0.
//
//    Input, int FACTOR_ORDER, the order of the factor.
//
//    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values for
//    factor FACTOR_INDEX.
//
//    Input, int FACTOR_NUM, the number of factors.
//
//    Input, int POINT_NUM, the number of elements in the direct product.
//
//    Input/output, double W[POINT_NUM], the elements of the
//    direct product, which are built up gradually.  
//
//  Local Parameters:
//
//    Local, integer START, the first location of a block of values to set.
//
//    Local, integer CONTIG, the number of consecutive values to set.
//
//    Local, integer SKIP, the distance from the current value of START
//    to the next location of a block of values to set.
//
//    Local, integer REP, the number of blocks of values to set.
//
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 1.0;
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( j = 0; j < factor_order; j++ )
  {
    start = 0 + j * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( i = start; i < start + contig; i++ )
      {
        w[i] = w[i] * factor_value[j];
      }
      start = start + skip;
    }
  }

  contig = contig * factor_order;

  return;
}
//****************************************************************************80

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_SOME prints "some" of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, integer I_LO, I_HI, the first and last indices to print.
//    The routine expects 1 <= I_LO <= I_HI <= N.
//
//    Input, string TITLE, an optional title.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    cout << "  " << setw(8)  << i       
         << "  " << setw(14) << a[i-1]  << "\n";
  }

  return;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n ) 
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

void sparse_grid_herm ( int dim_num, int level_max, int point_num, 
  double grid_weight[], double grid_point[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_HERM computes a sparse grid of Gauss-Hermite points.
//
//  Discussion:
//
//    The quadrature rule is associated with a sparse grid derived from 
//    a Smolyak construction using a 1D Gauss-Hermite quadrature rule.  
// 
//    The user specifies: 
//    * the spatial dimension of the quadrature region, 
//    * the level that defines the Smolyak grid. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, controls the size of the final sparse grid.
//
//    Input, int POINT_NUM, the number of points in the grid, as determined
//    by SPARSE_GRID_HERM_SIZE.
//
//    Output, double GRID_WEIGHT[POINT_NUM], the weights.
//
//    Output, double GRID_POINT[DIM_NUM*POINT_NUM], the points.
//
{
  int coeff;
  int dim;
  int *grid_base2;
  int *grid_index2;
  int *grid_level;
  double *grid_point_temp;
  double *grid_weight2;
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int *order_1d;
  int order_nd;
  int order_max;
  int point;
  int point_num2;
  int point2;
  int point3;
  int t;

  for ( point = 0; point < point_num; point++ )
  {
    grid_weight[point] = 0.0;
  }
//
//  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
//
  point_num2 = 0;

  level_min = i4_max ( 0, level_max + 1 - dim_num );
  
  grid_base2 = new int[dim_num];
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
//  that adds up to LEVEL.
//
    more = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );
//
//  Transform each 1D level to a corresponding 1D order.
//  The relationship is the same as for other OPEN rules.
//  The GL rule differs from the other OPEN rules only in the nesting behavior.
//
      level_to_order_open ( dim_num, level_1d, order_1d );
      
      for ( dim = 0; dim < dim_num; dim++ )
      {
        grid_base2[dim] = ( order_1d[dim] - 1 ) / 2;
      }
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = i4vec_product ( dim_num, order_1d );
//
//  Compute the weights for this product grid.
//
      grid_weight2 = product_weight_herm ( dim_num, order_1d, order_nd );
//
//  Now determine the coefficient of the weight.
//
      coeff = i4_power ( -1, level_max - level ) 
        * choose ( dim_num - 1, level_max - level );
//
//  The inner (hidden) loop generates all points corresponding to given grid.
//  The grid indices will be between -M to +M, where 2*M + 1 = ORDER_1D(DIM).
//
      grid_index2 = multigrid_index_z ( dim_num, order_1d, order_nd );
//
//  Determine the first level of appearance of each of the points.
//  This allows us to flag certain points as being repeats of points
//  generated on a grid of lower level.  
//
//  This is SLIGHTLY tricky.
//
      grid_level = index_level_herm ( level, level_max, dim_num, order_nd, 
        grid_index2, grid_base2 );
//
//  Only keep those points which first appear on this level.
//
      for ( point = 0; point < order_nd; point++ )
      {
//
//  Either a "new" point (increase count, create point, create weight)
//
        if ( grid_level[point] == level )
        {
          hermite_abscissa ( dim_num, 1, grid_index2+point*dim_num, 
            grid_base2, grid_point+point_num2*dim_num );

          grid_weight[point_num2] = ( double ) ( coeff ) * grid_weight2[point];
	  
	  point_num2 = point_num2 + 1;
	}
//
//  or an already existing point (create point temporarily, find match,
//  add weight to matched point's weight).
//
        else
        {
	  grid_point_temp = new double[dim_num];
	  
          hermite_abscissa ( dim_num, 1, grid_index2+point*dim_num, 
            grid_base2, grid_point_temp );
                     
	  for ( point2 = 0; point2 < point_num2; point2++ )
	  {
	    point3 = point2;
	    for ( dim = 0; dim < dim_num; dim++ )
	    {
              if ( grid_point[dim+point2*dim_num] != grid_point_temp[dim] )
	      {
	        point3 = -1;
	        break;
              }
	    }
            if ( point3 == point2 )
	    {
	      break;
	    }
          }
          
          if ( point3 == -1 )
	  {
            cout << "\n";
            cout << "SPARSE_GRID_HERM - Fatal error!\n";
            cout << "  Could not match point.\n";
            exit ( 1 );
          }
          
          grid_weight[point3] = grid_weight[point3] + 
            ( double ) ( coeff ) * grid_weight2[point];        
        }
      }

      delete [] grid_index2;
      delete [] grid_level;
      delete [] grid_weight2;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] grid_base2;
  delete [] level_1d;
  delete [] order_1d;
  
  return;
}
//****************************************************************************80

void sparse_grid_herm_index ( int dim_num, int level_max, int point_num, 
  int grid_index [], int grid_base[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_HERM_INDEX indexes points in a Gauss-Hermite sparse grid.
//
//  Discussion:
//
//    The sparse grid is assumed to be formed from 1D Gauss-Hermite rules
//    of ODD order, which have the property that only the central abscissa,
//    X = 0.0, is "nested".
//
//    The necessary dimensions of GRID_INDEX can be determined by 
//    calling SPARSE_GRID_HERM_SIZE first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Input, int POINT_NUM, the total number of points in the grids.
//
//    Output, int GRID_INDEX[DIM_NUM*POINT_NUM], a list of 
//    point indices, representing a subset of the product grid of level 
//    LEVEL_MAX, representing (exactly once) each point that will show up in a
//    sparse grid of level LEVEL_MAX.
//
//    Output, int GRID_BASE[DIM_NUM*POINT_NUM], a list of 
//    the orders of the Gauss-Hermite rules associated with each point 
//    and dimension.
//
{
  int dim;
  int factor;
  int *grid_base2;
  int *grid_index2;
  int *grid_level;
  int h;
  int j;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int *order_1d;
  int order_nd;
  int point;
  int point_num2;
  int t;
//
//  The outer loop generates LEVELs from LEVEL_MIN to LEVEL_MAX.
//
  point_num2 = 0;

  level_min = i4_max ( 0, level_max + 1 - dim_num );
  
  grid_base2 = new int[dim_num];
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates the next partition LEVEL_1D(1:DIM_NUM)
//  that adds up to LEVEL.
//
    more = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );
//
//  Transform each 1D level to a corresponding 1D order.
//
      level_to_order_open ( dim_num, level_1d, order_1d );
      for ( dim = 0; dim < dim_num; dim++ )
      {
        grid_base2[dim] = ( order_1d[dim] - 1 ) / 2;
      }
//
//  The product of the 1D orders gives us the number of points in this grid.
//
      order_nd = i4vec_product ( dim_num, order_1d );
//
//  The inner (hidden) loop generates all points corresponding to given grid.
//
      grid_index2 = multigrid_index_z ( dim_num, order_1d, order_nd );
//
//  Determine the first level of appearance of each of the points.
//  This allows us to flag certain points as being repeats of points
//  generated on a grid of lower level.  
//
//  This is SLIGHTLY tricky.
//
      grid_level = index_level_herm ( level, level_max, dim_num, order_nd, 
        grid_index2, grid_base2 );
//
//  Only keep those points which first appear on this level.
//
      for ( point = 0; point < order_nd; point++ )
      {
        if ( grid_level[point] == level )
        {
          for ( dim = 0; dim < dim_num; dim++ )
          {
            grid_index[dim+point_num2*dim_num] =
              grid_index2[dim+point*dim_num];
	    grid_base[dim+point_num2*dim_num] = grid_base2[dim];
          }
          point_num2 = point_num2 + 1;
        }
      }

      delete [] grid_index2;
      delete [] grid_level;

      if ( !more )
      {
        break;
      }
    }
  }
  delete [] grid_base2;
  delete [] level_1d;
  delete [] order_1d;

  return;
}
//****************************************************************************80

int sparse_grid_herm_size ( int dim_num, int level_max )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_GRID_HERM_SIZE sizes a sparse grid of Gauss-Hermite points.
//
//  Discussion:
//
//    The grid is defined as the sum of the product rules whose LEVEL
//    satisfies:
//
//      LEVEL_MIN <= LEVEL <= LEVEL_MAX.
//
//    where LEVEL_MAX is user specified, and 
//
//      LEVEL_MIN = max ( 0, LEVEL_MAX + 1 - DIM_NUM ).
//
//    The grids are only very weakly nested, since Gauss-Hermite rules
//    only have the origin in common.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEVEL_MAX, the maximum value of LEVEL.
//
//    Output, int SPARSE_GRID_HERM_SIZE, the number of points in the grid.
//
{
  int dim;
  int h;
  int level;
  int *level_1d;
  int level_min;
  bool more;
  int num;
  int *order_1d;
  int order_nd;
  int point_num;
  int t;
//
//  Special case.
//
  if ( level_max == 0 )
  {
    point_num = 1;
    return point_num;
  }
//
//  The outer loop generates LEVELs from 0 to LEVEL_MAX.
//
  point_num = 0;

  level_min = i4_max ( 0, level_max + 1 - dim_num );
  
  level_1d = new int[dim_num];
  order_1d = new int[dim_num];

  for ( level = level_min; level <= level_max; level++ )
  {
//
//  The middle loop generates the next partition that adds up to LEVEL.
//
    more = false;
    h = 0;
    t = 0;

    for ( ; ; )
    {
      comp_next ( level, dim_num, level_1d, &more, &h, &t );
//
//  Transform each 1D level to a corresponding 1D order.
//
      level_to_order_open ( dim_num, level_1d, order_1d );

      for ( dim = 0; dim < dim_num; dim++ )
      {
//
//  If we can reduce the level in this dimension by 1 and
//  still not go below LEVEL_MIN.
//
        if ( level_min < level && 1 < order_1d[dim] )
        {
          order_1d[dim] = order_1d[dim] - 1;
        }
      }
      point_num = point_num + i4vec_product ( dim_num, order_1d );

      if ( !more )
      {
        break;
      }
    }
  }

  delete [] level_1d;
  delete [] order_1d;

  return point_num;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void vec_colex_next2 ( int dim_num, int base[], int a[], bool *more )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_COLEX_NEXT2 generates vectors in colex order.
//
//  Discussion:
//
//    The vectors are produced in colexical order, starting with
//
//    (0,        0,        ...,0),
//    (1,        0,        ...,0),
//     ...
//    (BASE(1)-1,0,        ...,0)
//
//    (0,        1,        ...,0)
//    (1,        1,        ...,0)
//    ...
//    (BASE(1)-1,1,        ...,0)
//
//    (0,        2,        ...,0)
//    (1,        2,        ...,0)
//    ...
//    (BASE(1)-1,BASE(2)-1,...,BASE(DIM_NUM)-1).
//
//  Examples:
//
//    DIM_NUM = 2,
//    BASE = { 3, 3 }
//
//    0   0
//    1   0
//    2   0
//    0   1
//    1   1
//    2   1
//    0   2
//    1   2
//    2   2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int BASE[DIM_NUM], the bases to be used in each dimension.
//    In dimension I, entries will range from 0 to BASE[I]-1.
//
//    Output, int A[DIM_NUM], the next vector.
//
//    Input/output, bool *MORE.  Set this variable false before
//    the first call.  On return, MORE is TRUE if another vector has
//    been computed.  If MORE is returned FALSE, ignore the output 
//    vector and stop calling the routine.
//
{
  int i;

  if ( !( *more ) )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = 0;
    }
    *more = true;
  }
  else
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i] = a[i] + 1;

      if ( a[i] < base[i] )
      {
        return;
      }
      a[i] = 0;
    }
    *more = false;
  }

  return;
}
