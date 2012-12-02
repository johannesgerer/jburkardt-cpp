# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <string>

using namespace std;

int main ( );
void initialState ( char vis_table[], char real_table[], int n );
void printTable ( char vis_table[], int n );
void getneighbors ( int neighbors[], int pos, int n );
bool isRealNeighbor ( int index, int pos, int n );
void placeMines ( char real_table[], int n, int mines );
void calculateNumbers ( char real_table[], int n );
char countMines ( char real_table[], int neighbors[], int pos, int n );
void openSafe ( char real_table[], char vis_table[], int n, int pos, int &counter );
void play ( char real_table[], char vis_table[], int n, int mines );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main function for the MINESWEEPER program.
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    Detelina Stoyanova
//
//  Local Parameters:
//
//    Local, int MINES, the number of mines to use.
//
//    Local, int N, the number of rows and columns in the grid.
//
//    Local, char REAL_TABLE[N*N], contains a symbolic map with
//    '*' for mines, 
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Local, char VIS_TABLE[N*N], contains a symbolic map of facts the
//    user has uncovered:
//    '-', a square which has not been explored;
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
{
  int n = 8;
  int mines = 10;
  char real_table[n*n];
  char vis_table[n*n];
//
//  Initialize the random number generator.
//
  srand ( time ( NULL) );
//
//  Initialize both arrays with '-'.
//
  initialState ( vis_table, real_table, n );
//
//  In the real grid, replace N occurrences of '-' by '*'.
//  This is where the mines are.
//
  placeMines ( real_table, n, mines );
//
//  Count the mines adjacent to each square.
//
  calculateNumbers ( real_table, n );
//
//  Start the game.
//
  play ( real_table, vis_table, n, mines );

  return 0;
}
//****************************************************************************80

void initialState ( char vis_table[], char real_table[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    INITIALSTATE initializes all the elements of the two arrays to -'s
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    Detelina Stoyanova
//
//  Parameters:
//
//    Output, char VIS_TABLE[N*N], contains a symbolic map of facts the
//    user has uncovered:
//    '-', a square which has not been explored;
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Output, char REAL_TABLE[N*N], contains a symbolic map with
//    '*' for mines, 
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Input, int N, the number of rows and columns in the grid.
//
{
  for ( int i = 0; i < n * n; i++ )
  {
    vis_table[i] = '-';
    real_table[i] = '-';
  }
  return;
}
//****************************************************************************80

void printTable ( char vis_table[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRINTTABLE prints the visible table.
//
//  Discussion:
//
//    This function is called after each move by the player, to show the
//    current state of the game.
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    Detelina Stoyanova
//
//  Parameters:
//
//    Input, char VIS_TABLE[N*N], contains a symbolic map of facts the
//    user has uncovered:
//    '-', a square which has not been explored;
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Input, int N, the number of rows and columns in the grid.
//
{
  cout << "\n";
  cout << "  ";
  for ( int i = 1; i <= n; i++ )
  {
    cout << " " << i;
  } 
  cout << endl;
  for ( int i = 0; i < n; i++ )
  {
    cout << setw(2) << i+1 << " ";
    for ( int j = 0; j < n; j++ )
    {
      cout << vis_table[n*i + j] << " ";
    }
    cout << endl;
  }
  return;
}
//****************************************************************************80

void placeMines ( char real_table[], int n, int mines )

//****************************************************************************80
//
//  Purpose:
//
//    PLACEMINES places MINES mines randomly in the table.
//
//  Discussion:
//
//    On input, the N by N table contains only '-' values.  On output,
//    MINES entries of the table, chosen at random, now contain the value '*'.
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    Detelina Stoyanova
//
//  Parameters:
//
//    Input/output, char REAL_TABLE[N*N], contains a symbolic map with
//    '*' for mines, 
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Input, int N, the number of rows and columns in the grid.
//
//    Input, int MINES, the number of mines to use.
//
{
  int mine;
  int num_mines = 0;

  while ( num_mines < mines )
  {
    mine = rand ( ) % ( n * n );

    if ( real_table[mine] != '*' )
    {
      real_table[mine] = '*';
      num_mines++;
    }
  }
  return;
}
//****************************************************************************80

void calculateNumbers ( char real_table[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    CALCULATENUMBERS calculates the number of mines adjacent to each square.
//
//  Discussion:
//
//    This function is used to update the REAL_TABLE array once the position
//    of the mines has been randomly assigned.
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    Detelina Stoyanova
//
//  Parameters:
//
//    Input/output, char REAL_TABLE[N*N], contains a symbolic map with
//    '*' for mines, 
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Input, int N, the number of rows and columns in the grid.
//
{
  int counter;
  int neighbors[8]; 
//
//  For each square that is not mined, look at the eight neighbors,
//  count how many of them are mines, and store that value in the map.
//
  for ( int i = 0; i < n*n; i++ )
  {
    if ( real_table[i] != '*' )
    { 
      getneighbors ( neighbors, i, n );

      real_table[i] = countMines ( real_table, neighbors, i, n );
    }
  }
  return;
}
//****************************************************************************80

void getneighbors ( int neighbors[], int pos, int n )

//****************************************************************************80
//
//  Purpose:
//
//    GETNEIGHBORS returns the indices of the 8 neighbors of a square.
//
//  Discussion:
//
//    This function stores the position in the real_table of the 8 possible 
//    neighbors for each square (at position pos) I am using a one dimmensional 
//    so I have to take the number of squares per line (which is n) into account.
//
//    Here is a diagram of the indices of the neighbors of the square at POS:
//
//      ... |...    |...    |...    | ...
//         -+-------+-------+-------+-
//      ... |POS-N-1|POS-N  |POS-N+1| ...
//         -+-------+-------+-------+-
//      ... |POS-1  |POS    |POS+1  | ...
//         -+-------+-------+-------+-
//      ... |POS+N-1|POS+N  |POS+N+1| ...
//         -+-------+-------+-------+-
//      ... |...    |...    |...    | ...
//
//    If POS is near the top, bottom, left or right of the grid, then some of
//    these neighbors don't actually exist.  But we leave that problem to the
//    ISREALNEIGHBOR function.
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    Detelina Stoyanova
//
//  Parameters:
//
//    Output, int NEIGHBORS[8], the indices of the neighbors in the
//    northwest, north, northeast, west, east, southwest, south and southeast.
//
//    Input, int POS, the index of the square being considered.
//
//    Input, int N, the number of rows and columns in the grid.
//
{
  neighbors[0] = ( pos - n ) - 1;
  neighbors[1] =   pos - n ;
  neighbors[2] = ( pos - n ) + 1;

  neighbors[3] =   pos       - 1;
  neighbors[4] =   pos       + 1;

  neighbors[5] = ( pos + n ) - 1;
  neighbors[6] = ( pos + n );
  neighbors[7] = ( pos + n ) + 1;

  return;
}
//****************************************************************************80

bool isRealNeighbor ( int index, int pos, int n )

//****************************************************************************80
//
//  Purpose:
//
//    ISREALNEIGHBOR checks which neighbors are "real neighbors". 
//
//  Discussion:
//
//    On an N by N grid, the square at location POS may have up to eight
//    neighbors, in the northwest, north, northeast, west, east, southwest,
//    south and southeast directions.
//
//    However, some of these neighbors will not actually exist if
//    * POS is in the first row of the grid;
//    * POS is in the last row of the grid;
//    * POS is in the first column of the grid;
//    * POS is in the last column of the grid.
//
//    This function checks whether a neighbor in the INDEX direction
//    does not exist because of the placement of POS.
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    Detelina Stoyanova
//
//  Parameters:
//
//    Input, int INDEX, indicates the neighbor to be checked.
//
//    Input, int POS, the index of the square.
//
//    Input, int N, the number of rows and columns in the grid.
//
{
//
//  Northwest, north, northeast neighbors don't exist 
//  if POS < N.
//
  if ( pos < n )
  {
    if ( index == 0 || index == 1 || index == 2 )
    {
      return false;
    }
  }
//
//  Southwest, south, and southeast neighbors don't exist 
//  if N * ( N - 1 ) <= POS.
//
  if ( n * ( n - 1 ) <= pos )
  {
    if ( index == 5 || index == 6 || index == 7 )
    {
      return false;
    }
  }
//
//  Northwest, west, and southwest neighbors don't exist 
//  if POS is divisible by N.
//
  if ( pos % n == 0 )
  {
    if ( index == 0 || index == 3 || index == 5 )
    {
      return false;
    }
  }
//
//  Northeast, east, and southeast neighbors don't exist 
//  if POS + 1 is divisible by N.
//
  if ( ( pos + 1 ) % n == 0 )
  {
    if ( index == 2 || index == 4 || index == 7 )
    {
      return false;
    }
  }
//
//  Apparently, the neighbor in the INDEX direction is legal.
//
  return true;
}
//****************************************************************************80

char countMines ( char real_table[], int neighbors[], int pos, int n )

//****************************************************************************80
//
//  Purpose:
//
//    COUNTMINES computes the number of mines that are neighbors to a square.
//
//  Discussion:
//
//    This function takes a square at position POS and its 8 neighbors and
//    counts how many of the neighbors are mines.  
//
//    It returns the character ' ' if there are no neighboring mines,
//    otherwise a digit between 1 to 8 corresponding to the number of mines.
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    Detelina Stoyanova
//
//  Parameters:
//
//    Input/output, char REAL_TABLE[N*N], contains a symbolic map with
//    '*' for mines, 
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Input, int NEIGHBORS[8], the indices of the neighbors in the
//    northwest, north, northeast, west, east, southwest, south and southeast.
//
//    Input, int POS, the index of the square being examined.
//
//    Input, int N, the number of rows and columns in the grid.
//
{
  int mines = 0;
  char value;
//
//  Count the neighbor mines.
//
  for ( int i = 0; i < 8; i++ )
  { 
    if ( isRealNeighbor( i, pos, n ) )
    {
      if ( real_table[neighbors[i]] == '*' )
      {
        mines++;
      }
    }
  }
//
//  Form the character that symbolizes the neighbor mine count.
//
  if ( mines == 0 )
  {
   value = ' ';
  }
  else
  {
    value = ( char ) ( ( ( int ) '0' ) + mines );
  }
  return value;
}
//****************************************************************************80

void play ( char real_table[], char vis_table[], int n, int mines )

//****************************************************************************80
//
//  Purpose:
//
//    PLAY controls the game and allows the player to make moves.
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    Detelina Stoyanova
//
//  Parameters:
//
//    Input, char REAL_TABLE[N*N], contains a symbolic map with
//    '*' for mines, 
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Input/output, char VIS_TABLE[N*N], contains a symbolic map of facts the
//    user has uncovered:
//    '-', a square which has not been explored;
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Input, int N, the number of rows and columns in the grid.
//
//    Input, int MINES, the number of mines to use.
//
//  Local parameters:
//
//    Local, int COL, the column selected by the user.  This is a 1-based index.
//
//    Local, int COUNTER, counts the number of safe squares which the player
//    has cleared so far.  The game is over if the player can clear all the
//    safe squares.
//
//    Local, char CURR_SYMBOL, the symbol of the square which is currently
//    being cleared.  It could be blank, or a number or a '*'.
//
//    Local, int POS, the position in the 1D array corresponding to the
//    user's chosen ROW and COL.
//
//    Local, int ROW, the row selected by the user.  This is a 1-based index.
//
//    Local, string S, used to read and discard illegal input.
//
//    Local, int SAFE, the total number of safe squares.
{
  int col;
  int counter = 0;
  char curr_symbol;
  int pos;
  int row;
  string s;
  int safe = n*n - mines;
//
//  For each iteration of this loop, the player chooses a square, and the program
//  updates the game correspondingly.
//
  do
  {
//
//  Display the table of data the user has uncovered.
//
    printTable ( vis_table, n );   
//
//  Ask the player to specify the next square to clear.
//
    cout << "Enter row and column: " << endl;
    cin >> row 
        >> col;
//
//  Deal with ^D or end of file.
//
    if ( cin.eof ( ) )
    {
      cerr << "\n";
      cerr << "PLAY:\n";
      cerr << "  Input terminated with end-of-file.\n";
      cerr << "  Game terminated.\n";
      exit ( 1 );
    }
//
//  Deal with illegal input, such as "1,2" or "Quit".
//  You have to clear the error, then read the line as a string so it's gone.
//
    if ( cin.fail ( ) )
    {
      cin.clear ( );
      cin >> s;

      if ( s[0] == 'Q' )
      {
        cout << "\n";
        cout << "QUIT the game.\n";
        cout << "\n";
        cout << "Here is the mine map:\n";
        cout << "\n";
        printTable ( real_table, n );
        return;
      }

      if ( s[0] == 'M' )
      {
        cout << "\n";
        cout << "Enter row and column for mine flag: " << endl;
        cin >> row >> col;
        pos = n * ( row - 1 ) + ( col - 1 );
        if ( vis_table[pos] == '-' )
        {
          vis_table[pos] = '?';
        }
        cout << "\n";
        printTable ( vis_table, n );
        continue;
      }

      if ( s[0] == 'U' )
      {
        cout << "\n";
        cout << "Enter row and column for mine flag: " << endl;
        cin >> row >> col;
        pos = n * ( row - 1 ) + ( col - 1 );
        if ( vis_table[pos] == '?' )
        {
          vis_table[pos] = '-';
        }
        cout << "\n";
        printTable ( vis_table, n );
        continue;
      }

      cerr << "\n";
      cerr << "PLAY - Warning!\n";
      cerr << "  That input was illegal.\n";
      cerr << "  Please enter a row number and a column number.\n";

      continue;
    }

    if ( row < 1 || n < row || col < 1 || n < col )
    {
      cerr << "\n";
      cerr << "Your row or column index was illegal.\n";
      cerr << "Legal values are between 1 and " << n << ".\n";
      continue;
    }
//
//  Calculate the position of the square in the 1D array.
//
    pos = n * ( row - 1 ) + ( col - 1 );
//
//  If this square was not cleared already...
//
    if ( vis_table[pos] == '-' )
    {
//
//  ...then update this position in VIS_TABLE to equal the character which 
//  REAL_TABLE holds.
//
      curr_symbol = real_table[pos];
      vis_table[pos] = curr_symbol;
//
//  If the square which the player cleared has no mines around it
//  call a function which will automatically clear all its neighbors.
//  The function also increments the counter for the safe square which were open.
//
      if ( curr_symbol == ' ' )
      {
        openSafe ( real_table, vis_table, n, pos, counter );
      }
//
//  If the square was not a mine, we cleared another square, 
//  so increment the counter.
//
      if ( curr_symbol != '*' )
      {       
        counter++;
      }
    }

  } while ( counter < safe && curr_symbol != '*' ); 
//
//  If the number of squares cleared by the player is the same as 
//  the number of safe squares, then the player won the game.
//
  if ( counter == safe )
  {
    printTable ( real_table, n );
    cout << "\n";
    cout << "You Win!" << endl;
  }
//
//  The player loses by choosing a mine.
//  The position of all mines is shown.
//
  if ( curr_symbol == '*' )
  {
    for ( int i = 0; i < n * n; i++ )
    {
      if ( real_table[i] == '*' )
      {
        vis_table[i] = '*';
      }
    }
    printTable ( vis_table, n );
    cout << "\n";
    cout << "You Lose!" << endl;
  }
  return;
}
//****************************************************************************80

void openSafe ( char real_table[], char vis_table[], int n, int pos, 
  int &counter )

//****************************************************************************80
//
//  Purpose:
//
//    OPENSAFE is called whenever a square with no neighboring mines is cleared.
//
//  Discussion:
//
//    If there are no mines around the square at position POS then it is safe 
//    to clear all the squares around it.  It is possible that some of the 
//    neighbors also don't have mines around them so their neighbors have to 
//    get cleared as well.  This is why this function uses recursion (i.e it 
//    calls itself when a neighbors has 0 mines around it).  
//
//  Modified:
//
//    24 July 2011
//
//  Author:
//
//    Detelina Stoyanova
//
//  Parameters:
//
//    Input, char REAL_TABLE[N*N], contains a symbolic map with
//    '*' for mines, 
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Local, char VIS_TABLE[N*N], contains a symbolic map of facts the
//    user has uncovered:
//    '-', a square which has not been explored;
//    ' ' for unmined squares with no mine neighbors,
//    a digit between 1 and 8 for unmined squares with some mine neighbors.
//
//    Input, int N, the number of rows and columns in the grid.
//
//    Input, int POS, the index of the square being considered.
//
{
  int neighbors[8];
//
//  Index the 8 neighbors of this square.
//
  getneighbors ( neighbors, pos, n );
//
//  Consider each neighbor square.
//
  for ( int i = 0; i < 8; i++ )
  { 
//
//  If there really is a neighbor square in that direction, ...
//
    if ( isRealNeighbor ( i, pos, n ) )
    {
//
//  ... and the neighbor square is not a mine, and has not been explored yet ...
//
      if ( vis_table[neighbors[i]] == '-' )
      {
//
//  ...THEN make the neighbor square visible and add one more to the
//  count of squares cleared.
//
        vis_table[neighbors[i]] = real_table[neighbors[i]];
        counter++;
//
//  If the neighbor square has no mine neighbors, then we have to
//  check the neighbors of this neighbor, now, too!
//
        if ( vis_table[neighbors[i]] == ' ' )
        {
          openSafe ( real_table, vis_table, n, neighbors[i], counter );
        }
      }
    } 
  }
  return;
}

