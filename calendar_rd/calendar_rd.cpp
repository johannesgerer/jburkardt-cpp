// The following C++ code is translated from the Lisp code
// in ``Calendrical Calculations'' by Nachum Dershowitz and
// Edward M. Reingold, Software---Practice & Experience,
// vol. 20, no. 9 (September, 1990), pp. 899--928.

// This code is in the public domain, but any use of it
// should publicly acknowledge its source.

// Classes GregorianDate, JulianDate, IsoDate, IslamicDate,
// and HebrewDate

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <cmath>

using namespace std;

class IsoDate;

char *DayName[7] = 
{
  "Sunday", 
  "Monday", 
  "Tuesday", 
  "Wednesday",
  "Thursday", 
  "Friday", 
  "Saturday"
};
//
//  Absolute dates of the starts of the Hebrew, Islamic, and Julian calendars.
//
const int HebrewEpoch = -1373429;
const int IslamicEpoch = 227014;
const int JulianEpoch = -2;

// Absolute dates

// "Absolute date" means the number of days elapsed since the Gregorian date
// Sunday, December 31, 1 BC. (Since there was no year 0, the year following
// 1 BC is 1 AD.) Thus the Gregorian date January 1, 1 AD is absolute date
// number 1.

//****************************************************************************80

int XdayOnOrBefore ( int d, int x ) 

//****************************************************************************80
//
//  Purpose:
//
//    XdayOnOrBefore returns the date of the weekday before the given date.
//
//  Modified:
//
//    07 December 2008
//
//  Author:
//
//    Nachum Dershowitz, Edward Reingold
//
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz,
//    Calendrical Calculations I,
//    Software - Practice and Experience,
//    Volume 20, Number 9, September 1990, pages 899-928.
//
//  Parameters:
//
//    Input, int D, the absolute date.
//
//    Input, int X, the weekday (0 through 6, with 0 meaning Sunday).
//
//    Output, int XdayOnOrBefore, the date of the last weekday of 
//    type X on or before date D.
//
{
  return (d - ((d - x) % 7));
}
//****************************************************************************80

int LastDayOfGregorianMonth ( int month, int year ) 

//****************************************************************************80
//
//  Purpose:
//
//    LastDayOfGregorianMonth computes the last date of the month for the Gregorian calendar.
//
//  Modified:
//
//    07 December 2008
//
//  Author:
//
//    Nachum Dershowitz, Edward Reingold
//
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz,
//    Calendrical Calculations I,
//    Software - Practice and Experience,
//    Volume 20, Number 9, September 1990, pages 899-928.
//
//  Parameters:
//
{  
  switch (month) {
  case 2:
    if ((((year % 4) == 0) && ((year % 100) != 0))
        || ((year % 400) == 0))
      return 29;
    else
      return 28;
  case 4:
  case 6:
  case 9:
  case 11: return 30;
  default: return 31;
  }
}

class GregorianDate {
private:
  int year;   // 1...
  int month;  // 1 == January, ..., 12 == December
  int day;    // 1..LastDayOfGregorianMonth(month, year)
  
public:
  GregorianDate(int m, int d, int y) { month = m; day = d; year = y; }
  
  GregorianDate(int d) { // Computes the Gregorian date from the absolute date.
    
    // Search forward year by year from approximate year
    year = d/366;
    while (d >= GregorianDate(1,1,year+1))
      year++;
    // Search forward month by month from January
    month = 1;
    while (d > GregorianDate(month, LastDayOfGregorianMonth(month,year), year))
      month++;
    day = d - GregorianDate(month,1,year) + 1;
  }
  
  operator int() { // Computes the absolute date from the Gregorian date.
    int N = day;           // days this month
    for (int m = month - 1;  m > 0; m--) // days in prior months this year
      N = N + LastDayOfGregorianMonth(m, year);
    return
      (N                    // days this year
       + 365 * (year - 1)   // days in previous years ignoring leap days
       + (year - 1)/4       // Julian leap days before this year...
       - (year - 1)/100     // ...minus prior century years...
       + (year - 1)/400);   // ...plus prior years divisible by 400
  }
  
  int GetMonth() { return month; }
  int GetDay() { return day; }
  int GetYear() { return year; }
  
};

ostream& operator<<(ostream& c, GregorianDate d) {
  c << d.GetMonth() << "/" << d.GetDay() << "/" << d.GetYear();
  return c;
};

//****************************************************************************80

GregorianDate NthXday(int n, int x, int month, int year, int day = 0)

//****************************************************************************80
//
// The Gregorian date of nth x-day in month, year before/after optional day.
// x = 0 means Sunday, x = 1 means Monday, and so on.  If n<0, return the nth
// x-day before month day, year (inclusive).  If n>0, return the nth x-day
// after month day, year (inclusive).  If day is omitted or 0, it defaults
// to 1 if n>0, and month's last day otherwise.
{
  if (n > 0) {
    if (day == 0)
      day = 1;  // default for positive n
    return GregorianDate
      ((7 * (n - 1)) + XdayOnOrBefore(6 + GregorianDate(month, day, year), x));
  }
  else {
    if (day == 0)
      day = LastDayOfGregorianMonth(month, year);;  // default for negative n
    return GregorianDate
      ((7 * (n + 1)) + XdayOnOrBefore(GregorianDate(month, day, year), x));
  }
}
//****************************************************************************80

int LastDayOfJulianMonth(int month, int year) 

//****************************************************************************80
{
// Compute the last date of the month for the Julian calendar.
  switch (month) {
  case 2:
    if ((year % 4) == 0)
      return 29;
    else
      return 28;
  case 4:
  case 6:
  case 9:
  case 11: return 30;
  default: return 31;
  }
}

class JulianDate {
private:
  int year;   // 1...
  int month;  // 1 == January, ..., 12 == December
  int day;    // 1..LastDayOfJulianMonth(month, year)
  
public:
  JulianDate(int m, int d, int y) { month = m; day = d; year = y; }
  
  JulianDate(int d) { // Computes the Julian date from the absolute date.
    // Search forward year by year from approximate year
    year = (d + JulianEpoch)/366;
    while (d >= JulianDate(1,1,year+1))
      year++;
    // Search forward month by month from January
    month = 1;
    while (d > JulianDate(month, LastDayOfJulianMonth(month,year), year))
      month++;
    day = d - JulianDate(month,1,year) + 1;
  }
  
  operator int() { // Computes the absolute date from the Julian date.
    
    int N = day;                         // days this month
    for (int m = month - 1;  m > 0; m--) // days in prior months this year
      N = N + LastDayOfJulianMonth(m, year);
    return
      (N                     // days this year
       + 365 * (year - 1)    // days in previous years ignoring leap days
       + (year - 1)/4        // leap days before this year...
       + JulianEpoch);       // days elapsed before absolute date 1
  }
  
  int GetMonth() { return month; }
  int GetDay() { return day; }
  int GetYear() { return year; }
  int absolute();
  
};

ostream& operator<<(ostream& c, JulianDate d) {
  c << d.GetMonth() << "/" << d.GetDay() << "/" << d.GetYear();
  return c;
};


// ISO dates

class IsoDate {
private:
  int year;  // 1...
  int week;  // 1..52 or 53
  int day;   // 1..7
  
public:
  IsoDate(int w, int d, int y) { week = w; day = d; year = y; }
  
  IsoDate(int d) { // Computes the ISO date from the absolute date.
    year = GregorianDate(d - 3).GetYear();
    if (d >= IsoDate(1,1,year+1))
      year++;
    if ((d % 7) == 0)
      day = 7;      // Sunday
    else
      day = d % 7;  // Monday..Saturday
    week = 1 + (d - IsoDate(1,1,year)) / 7;
  }
  
  operator int() { // Computes the absolute date from the ISO date.
    return
      XdayOnOrBefore(GregorianDate(1,4,year),1) // days in prior years
      + 7 * (week - 1)                          // days in prior weeks this year
      + (day - 1);                              // prior days this week
  }
  
  int GetWeek() { return week; }
  int GetDay() { return day; }
  int GetYear() { return year; }
  
};

ostream& operator<<(ostream& c, IsoDate d) {
  c << d.GetWeek() << "/" << d.GetDay() << "/" << d.GetYear();
  return c;
}
//****************************************************************************80

int IslamicLeapYear(int year) 

//****************************************************************************80
{
// True if year is an Islamic leap year
  
  if ((((11 * year) + 14) % 30) < 11)
    return 1;
  else
    return 0;
}
//****************************************************************************80

int LastDayOfIslamicMonth(int month, int year) 

//****************************************************************************80
{
// Last day in month during year on the Islamic calendar.
  
  if (((month % 2) == 1) || ((month == 12) && IslamicLeapYear(year)))
    return 30;
  else
    return 29;
}

class IslamicDate {
private:
  int year;   // 1...
  int month;  // 1..13 (12 in a common year)
  int day;    // 1..LastDayOfIslamicMonth(month,year)
  
public:
  IslamicDate(int m, int d, int y) { month = m; day = d; year = y; }
  
  IslamicDate(int d) { // Computes the Islamic date from the absolute date.
    if (d <= IslamicEpoch) { // Date is pre-Islamic
      month = 0;
      day = 0;
      year = 0;
    }
    else {
      // Search forward year by year from approximate year
      year = (d - IslamicEpoch) / 355;
      while (d >= IslamicDate(1,1,year+1))
        year++;
      // Search forward month by month from Muharram
      month = 1;
      while (d > IslamicDate(month, LastDayOfIslamicMonth(month,year), year))
        month++;
      day = d - IslamicDate(month,1,year) + 1;
    }
  }
  
  operator int() { // Computes the absolute date from the Islamic date.
    return (day                      // days so far this month
            + 29 * (month - 1)       // days so far...
            + month/2                //            ...this year
            + 354 * (year - 1)       // non-leap days in prior years
            + (3 + (11 * year)) / 30 // leap days in prior years
            + IslamicEpoch);                // days before start of calendar
  }
  
  int GetMonth() { return month; }
  int GetDay() { return day; }
  int GetYear() { return year; }
  
};

ostream& operator<<(ostream& c, IslamicDate d) {
  c << d.GetMonth() << "/" << d.GetDay() << "/" << d.GetYear();
  return c;
}
//****************************************************************************80

int HebrewLeapYear(int year) 

//****************************************************************************80
{
// True if year is an Hebrew leap year
  
  if ((((7 * year) + 1) % 19) < 7)
    return 1;
  else
    return 0;
}
//****************************************************************************80

int LastMonthOfHebrewYear(int year) 

//****************************************************************************80
{
// Last month of Hebrew year.
  
  if (HebrewLeapYear(year))
    return 13;
  else
    return 12;
}
//****************************************************************************80

int HebrewCalendarElapsedDays(int year) 

//****************************************************************************80
{
// Number of days elapsed from the Sunday prior to the start of the
// Hebrew calendar to the mean conjunction of Tishri of Hebrew year.
  
  int MonthsElapsed =
    (235 * ((year - 1) / 19))           // Months in complete cycles so far.
    + (12 * ((year - 1) % 19))          // Regular months in this cycle.
    + (7 * ((year - 1) % 19) + 1) / 19; // Leap months this cycle
  int PartsElapsed = 204 + 793 * (MonthsElapsed % 1080);
  int HoursElapsed =
    5 + 12 * MonthsElapsed + 793 * (MonthsElapsed  / 1080)
    + PartsElapsed / 1080;
  int ConjunctionDay = 1 + 29 * MonthsElapsed + HoursElapsed / 24;
  int ConjunctionParts = 1080 * (HoursElapsed % 24) + PartsElapsed % 1080;
  int AlternativeDay;
  if ((ConjunctionParts >= 19440)        // If new moon is at or after midday,
      || (((ConjunctionDay % 7) == 2)    // ...or is on a Tuesday...
          && (ConjunctionParts >= 9924)  // at 9 hours, 204 parts or later...
          && !(HebrewLeapYear(year)))   // ...of a common year,
      || (((ConjunctionDay % 7) == 1)    // ...or is on a Monday at...
          && (ConjunctionParts >= 16789) // 15 hours, 589 parts or later...
          && (HebrewLeapYear(year - 1))))// at the end of a leap year
    // Then postpone Rosh HaShanah one day
    AlternativeDay = ConjunctionDay + 1;
  else
    AlternativeDay = ConjunctionDay;
  if (((AlternativeDay % 7) == 0)// If Rosh HaShanah would occur on Sunday,
      || ((AlternativeDay % 7) == 3)     // or Wednesday,
      || ((AlternativeDay % 7) == 5))    // or Friday
    // Then postpone it one (more) day
    return (1+ AlternativeDay);
  else
    return AlternativeDay;
}
//****************************************************************************80

int DaysInHebrewYear(int year) 

//****************************************************************************80
{
// Number of days in Hebrew year.
  
  return ((HebrewCalendarElapsedDays(year + 1)) -
          (HebrewCalendarElapsedDays(year)));
}
//****************************************************************************80

int LongHeshvan(int year) 

//****************************************************************************80
{
// True if Heshvan is long in Hebrew year.
  
  if ((DaysInHebrewYear(year) % 10) == 5)
    return 1;
  else
    return 0;
}
//****************************************************************************80

int ShortKislev(int year) 

//****************************************************************************80
{
// True if Kislev is short in Hebrew year.
  
  if ((DaysInHebrewYear(year) % 10) == 3)
    return 1;
  else
    return 0;
}
//****************************************************************************80

int LastDayOfHebrewMonth(int month, int year) 

//****************************************************************************80
{
// Last day of month in Hebrew year.
  
  if ((month == 2)
      || (month == 4)
      || (month == 6)
      || ((month == 8) && !(LongHeshvan(year)))
      || ((month == 9) && ShortKislev(year))
      || (month == 10)
      || ((month == 12) && !(HebrewLeapYear(year)))
      || (month == 13))
    return 29;
  else
    return 30;
}

class HebrewDate {
private:
  int year;   // 1...
  int month;  // 1..LastMonthOfHebrewYear(year)
  int day;    // 1..LastDayOfHebrewMonth(month, year)
  
public:
  HebrewDate(int m, int d, int y) { month = m; day = d; year = y; }
  
  HebrewDate(int d) { // Computes the Hebrew date from the absolute date.
    year = (d + HebrewEpoch) / 366; // Approximation from below.
    // Search forward for year from the approximation.
    while (d >= HebrewDate(7,1,year + 1))
      year++;
    // Search forward for month from either Tishri or Nisan.
    if (d < HebrewDate(1, 1, year))
      month = 7;  //  Start at Tishri
    else
      month = 1;  //  Start at Nisan
    while (d > HebrewDate(month, (LastDayOfHebrewMonth(month,year)), year))
      month++;
    // Calculate the day by subtraction.
    day = d - HebrewDate(month, 1, year) + 1;
  }
  
  operator int() { // Computes the absolute date of Hebrew date.
    int DayInYear = day; // Days so far this month.
    if (month < 7) { // Before Tishri, so add days in prior months
                     // this year before and after Nisan.
      int m = 7;
      while (m <= (LastMonthOfHebrewYear(year))) {
        DayInYear = DayInYear + LastDayOfHebrewMonth(m, year);
        m++;
      };
      m = 1;
      while (m < month) {
        DayInYear = DayInYear + LastDayOfHebrewMonth(m, year);
        m++;
      }
    }
    else { // Add days in prior months this year
      int m = 7;
      while (m < month) {
        DayInYear = DayInYear + LastDayOfHebrewMonth(m, year);
        m++;
      }
    }
    return (DayInYear +
            (HebrewCalendarElapsedDays(year)// Days in prior years.
             + HebrewEpoch));         // Days elapsed before absolute date 1.
  }
  
  int GetMonth() { return month; }
  int GetDay() { return day; }
  int GetYear() { return year; }
  
};

ostream& operator<<(ostream& c, HebrewDate d) {
  c << d.GetMonth() << "/" << d.GetDay() << "/" << d.GetYear();
  return c;
}
//****************************************************************************80

int main ( ) 

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CALENDAR_RD.
//
//  Modified:
//
//    18 January 2009
//
//  Author:
//
//    Nachum Dershowitz, Edward Reingold
//
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz,
//    Calendrical Calculations I,
//    Software - Practice and Experience,
//    Volume 20, Number 9, September 1990, pages 899-928.
//
{
  int d;
  double jed;
  int m;
  int y;
  
  while ( 1 ) 
  {
    cout << "Enter year (>0): ";                        cin >> y;
    if ( y <= 0 ) 
    {
      break;
    }

    cout << "Enter month (1..12): ";                    cin >> m;
    cout << "Enter day (1.."
         << LastDayOfGregorianMonth(m, y) << "): "; cin >> d;
    
    GregorianDate g(m,d,y);
    int a = g;
    cout << g << " = " << a << " = " << DayName[g % 7] << "\n";
    
    g = a;
    a = g;
    cout << "    = Gregorian date " << g << " = absolute date " << a << "\n";
    
    JulianDate j(a);
    a = j;
    cout << "    = Julian date " << j << " = absolute date " << a << "\n";
    
    IsoDate i(a);
    a = i;
    cout << "    = ISO date " << i << " = absolute date " << a << "\n";
    
    HebrewDate h(a);
    a = h;
    cout << "    = Hebrew date " << h << " = absolute date " << a << "\n";
    
    IslamicDate I(a);
    a = I;
    cout << "    = Islamic date " << I << " = absolute date " << a << "\n";  

    jed = a + 1721424.5;
    cout << "    = JED " << setw(20) << setprecision(2) << fixed << jed << "\n";
  }

  return 0;
}
