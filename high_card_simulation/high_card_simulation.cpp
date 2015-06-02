# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "high_card_simulation.hpp"

//****************************************************************************80

double *high_card_probability ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    HIGH_CARD_PROBABILITY: winning probabilities for the high card game.
//
//  Discussion:
//
//    The high card game presents the player with a deck of cards, each
//    having an unknown value.  The player is allowed to go throught the
//    deck once, looking at the cards one at a time.  At any time, the player
//    may decide to take a particular card, winning that amount and stopping
//    the game.  If the player continues to the end, by default the last card
//    indicates the amount won.
//
//    An optimal strategy for selecting the highest card is as follows:
//    * look at, but do not select, the first k-1 cards;
//    * stop at the first card, from k to n, that is higher than the first 
//      k-1 cards.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of cards.  
//
//    Output, double P[N].  P[K] is the probability that a strategy 
//    that skips K cards will win, given that the deck has N cards.
//
{
  int i;
  int j;
  double *p;
  double t;

  p = new double[n];

  for ( i = 0; i < n; i++ )
  {
    t = 0.0;
    for ( j = i + 1; j < n; j++ )
    {
      t = t + 1.0 / ( double ) ( j );
    }
    p[i] = ( 1.0 + ( double ) ( i ) * t ) / ( double ) ( n );
  }

  return p;
}
//****************************************************************************80

int *high_card_shuffle ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    HIGH_CARD_SHUFFLE generates a sequence of numeric "cards" for a game.
//
//  Discussion:
//
//    In this game, you know that the deck contains N cards.  You win by 
//    choosing the highest card in the deck.  You don't know what this card
//    is, and you must choose your card by saying "stop" as, one by one,
//    the cards of the deck are exposed.  
//
//    A random guesser would get the high card with probability 1/N.
//
//    An intelligent guesser can do much better.
//
//    It is the goal of this program so "shuffle" a deck of cards suitable
//    for this game.  The problem is that we know the highest card in an
//    ordinary deck.  Let's replace the cards by integers.  Then if we know
//    in advance the range of the cards (say, they must lie between 1 and 
//    1,000), it may be true that we can guess the card that is the maximum.
//
//    However, this program produces a sequence of integer card values for
//    which no information can be gained from the values.  It does this
//    by regarding the card values as binary integers between 1 and 2^N - 1.
//    We can make a perfectly information-free sequence as follows:
//
//      Card 1 sets bit N-1 to 1.
//      Card 2 sets bit N-2 to 1, bit  N-1 randomly.
//      ...
//      Card I sets bit N-I to 1, bits N-1 down to N-I+1 randomly.
//      ...
//      Card N sets bit N-N to 1, bits N-1 down to 1 randomly.
//
//    The I-th card has equal probability to fall in any of the I intervals
//    defined by the I-1 previous cards.  So, knowing the previous cards tells
//    you absolutely nothing about where the next card will fall, and each
//    card is, at the moment you see it, as good a guess for the maximum as
//    the unseen cards.
//
//    For example, the command "high_card_shuffle(7)" might yield
//
//      64    96    80     8     4    82    29
//    or
//      64    32    16    24    12    58    73
//    or
//      64    96    48     8   100    26    93
//    or
//      64    96    16    56    76   122    71
//
//    in which the highest card is #2, #7, #5, or #6 respectively.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of cards.  N probably needs to 
//    be less than 32.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, int SEQUENCE[N], a set of N integer values 
//    that can be used as the cards in the high card guessing game.
//
{
  int c;
  int i;
  int j;
  int k;
  int *sequence;

  if ( 32 <= n )
  {
    cerr << "\n";
    cerr << "HIGH_CARD_SHUFFLE - Fatal error!\n";
    cerr << "  This program can only handle N < 32.\n";
    exit ( 1 );
  }

  sequence = new int[n];

  for ( i = 0; i < n; i++ )
  {
    c = i4_power ( 2, n - i - 1 );
    for ( j = 0; j < i; j++ )
    {
      k = i4_uniform_ab ( 0, 1, seed );
      c = c + k * i4_power ( 2, n - i + j );
    }
    sequence[i] = c;
  }

  return sequence;
}
//****************************************************************************80

double high_card_simulation ( int deck_size, int skip_num, int trial_num,
  int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    HIGH_CARD_SIMULATION simulates a game of choosing the highest card in a deck.
//
//  Discussion:
//
//    You are given a deck of DECK_SIZE cards.
//
//    Your goal is to select the high card.  For convenience, we can assume 
//    the cards are a permutation of the integers from 1 to DECK_SIZE, but in
//    fact the user mustn't see such values or else it's obvious which is the
//    largest card.
//
//    However, your choice is made under the following rules:  You may turn over
//    one card at a time.  When a card is turned over, you may declare that to be
//    your choice, or else turn over another card.  If you have not chosen a card
//    by the end, then your choice is the final card.
//
//    If you have no idea what to do, and simply decide in advance to pick
//    a card "at random", that is, for example, you decide to pick the 15th card
//    before having seen any cards, then your probability of winning is 
//    1/DECK_SIZE.
//
//    The question is, can you do better than that?
//
//    Your strategy is as follows: always look at the first SKIP_NUM cards 
//    without choosing them.  Then choose the very next card you encounter 
//    that is larger than the cards you skipped.
//
//    Using this program, you can easily see that skipping 5 cards is much better
//    than picking one at random, skipping 10 is even better, and so on.  
//    Of course, you can't skip too many cards, and in fact, the results seem 
//    to be best for somewhere around 30 to 35 cards skipped.  For problems 
//    like this, the optimal value is somewhere around 1 / e, where E is the 
//    base of the natural logarithm system.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DECK_SIZE, the number of cards in the deck.
//    2 <= DECK_SIZE.  Default value is 52;
//
//    Input, int SKIP_NUM, the number of initial cards you plan 
//    to examine but will NOT select.  If SKIP_NUM is 0, you don't look at any 
//    cards first.  0 <= SKIP_NUM < DECK_SIZE.
//
//    Input, int TRIAL_NUM, the number of times we will 
//    simulate this process.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double HIGH_CARD_SIMULATION, the estimated probability that 
//    your strategy of skipping SKIP_NUM cards and then selecting the next 
//    card that is bigger, will result in choosing the highest card.
//
{
  int card;
  int *cards;
  int choice;
  int correct;
  int i4_huge = 2147483647;
  double p;
  int skip_max;
  int trial;
  int true_max;
//
//  Check values.
//
  if ( deck_size < 2 )
  {
    cerr << "\n";
    cerr << "HIGH_CARD_SIMULATION - Fatal error!\n";
    cerr << "  DECK_SIZE must be at least 2.\n";
    cerr << "  Your value was %d\n", deck_size;
    exit ( 1 );
  }

  if ( skip_num < 0 )
  {
    skip_num = 0;
  }

  if ( deck_size <= skip_num )
  {
    cerr << "\n";
    cerr << "HIGH_CARD_SIMULATION - Fatal error!\n";
    cerr << "  SKIP_NUM must be less than DECK_SIZE.\n";
    cerr << "  Your DECK_SIZE = " << deck_size << "\n";
    cerr << "  Your SKIP_NUM = " << skip_num << "\n";
    exit ( 1 );
  }

  if ( trial_num < 1 )
  {
    cerr << "\n";
    cerr << "HIGH_CARD_SIMULATION - Fatal error!\n";
    cerr << "  TRIAL_NUM must be at least 1.\n";
    cerr << "  Your TRIAL_NUM was " << trial_num << "\n";
    exit ( 1 );
  }

  correct = 0;

  for ( trial = 1; trial <= trial_num; trial++ )
  {
    cards = perm_uniform_new ( deck_size, seed );

    if ( 1 <= skip_num )
    {
      skip_max = i4vec_max ( skip_num, cards );
    }
    else
    {
      skip_max = - i4_huge;
    }

    true_max = i4vec_max ( deck_size, cards );
//
//  In case you don't encounter a card larger than SKIP_MAX,
//  we'll assume you pick the last card in the deck, even though
//  you know it's a loser.
//
    choice = cards[deck_size-1];
//
//  Turn over the remaining cards in the deck, but stop
//  immediately when you find one bigger than SKIP_MAX.
//
    for ( card = skip_num; card < deck_size; card++ )
    {
      if ( skip_max < cards[card] )
      {
        choice = cards[card];
        break;
      }
    }
//
//  Record successful choices.
//
    if ( choice == true_max )
    {
      correct = correct + 1;
    }
    free ( cards );
  }
//
//  Estimate the probability.
//
  p = ( double ) ( correct ) / ( double ) ( trial_num );

  return p;
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
//    01 April 2004
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
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
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
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
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

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

int i4vec_max ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MAX returns the value of the maximum element in an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MAX, the value of the maximum element.  This
//    is set to 0 if N <= 0.
//
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < a[i] )
    {
      value = a[i];
    }
  }

  return value;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8) << i
         << ": " << setw(8) << a[i]  << "\n";
  }
  return;
}
/******************************************************************************/

int *perm_uniform_new ( int n, int &seed )

/******************************************************************************/
//
//  Purpose:
//
//    PERM_UNIFORM_NEW selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int PERM_UNIFORM_NEW[N], a permutation of
//    (BASE, BASE+1, ..., BASE+N-1).
//
{
  int i;
  int j;
  int k;
  int *p;

  p = new int[n];

  for ( i = 0; i < n; i++ )
  {
    p[i] = i;
  }

  for ( i = 0; i < n - 1; i++ )
  {
    j = i4_uniform_ab ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }

  return p;
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
