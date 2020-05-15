/***********************************************************
 *														   *
 *   Martin Boer											   *
 *   Biometris Wageningen-UR & NWO, grant no. 925-01-001    *
 *   July 12, 2002										   *
 *														   *
 *   last update : 25 aug 2003							   *
 ************************************************************/
#ifndef UTIL_HEADER
#define UTIL_HEADER

//#pragma warning(disable:4786)   // disable C4786 warning

#include <math.h>
#include <time.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "mblexcept.h"
#include "misc.h"
#include "convert.h"

namespace mbl
{

// const double ONE_PI			 = 3.141592653589793;
// const double TWO_PI			 = 6.283185307179586;
// const double INV_SQRT_TWO_PI = 0.3989422804014326779399461;
// const double LOG10E			 = 0.4342944819032518276511289;
//
// inline double abs(double x) { return fabs(x);}

// double mean(const std::vector<double>& y);
// double mean(const std::vector<double>& y, const std::vector<double>& w);
//
// double variance(const std::vector<double>& y);
// double variance(const std::vector<double>& y, const std::vector<double>& w);
// double covariance(const std::vector<double>& x, const std::vector<double>& y);
//
// double sum_of_sqr(const std::vector<double>& x);
// double sum_of_sqr(const std::vector<double>& x, const std::vector<double>& w);

/*
 template<class T1, class T2>
 double inner_product(const std::vector<T1>& a, const std::vector<T2>& b)
 {
 std::vector<T1>::const_iterator iter1;
 std::vector<T2>::const_iterator iter2;

 double sum = 0.0;
 for (iter1 = a.begin(), iter2 = b.begin();
 iter1 != a.end(); iter1++, iter2++)
 sum += (*iter1)*(*iter2);
 return sum;
 }
 */

// inline double norm_sqr(const std::vector<double>& x) { return sum_of_sqr(x); }
// inline double norm(const std::vector<double>& x) { return sqrt(norm_sqr(x)); }


// template<class T1> void read_columns(std::vector<T1>& x,
//                                      int dim,
//                                      std::istream& inp)
// {
//   std::string line;
//   for (int i=0;i<dim;i++)
//   {
//     getline(inp,line);
//     std::istringstream line_stream(line);
//     T1 val_x;
//     line_stream >> val_x;
//     x.push_back(val_x);
//   }
// }
//
// template<class T1,class T2> void read_columns(std::vector<T1>& x,
//                                               std::vector<T2>& y,
//                                               int dim,
//                                               std::istream& inp)
// {
//   std::string line;
//   for (int i=0;i<dim;i++)
//   {
//     getline(inp,line);
//     std::istringstream line_stream(line);
//     T1 val_x;
//     T2 val_y;
//     line_stream >> val_x >> val_y;
//     x.push_back(val_x);
//     y.push_back(val_y);
//   }
// }
//
// template<class T1,class T2,class T3> void read_columns(std::vector<T1>& x,
//                                                        std::vector<T2>& y,
//                                                        std::vector<T3>& z,
//                                                        int dim,
//                                                        std::istream& inp)
// {
//   std::string line;
//   for (int i=0;i<dim;i++)
//   {
//     getline(inp,line);
//     std::istringstream line_stream(line);
//     T1 val_x;
//     T2 val_y;
//     T3 val_z;
//     line_stream >> val_x >> val_y >> val_z;
//     x.push_back(val_x);
//     y.push_back(val_y);
//     z.push_back(val_z);
//   }
// }

// template<class InputIterator>
// bool all_equal(InputIterator first, InputIterator last)
// {
// 	InputIterator iter = first;
// 	++iter;
// 	for (; iter != last; ++iter)
// 		if (*iter != *first)
// 			return false;
// 	return true;
// }

// class Stopwatch
// {
// public:
// 	Stopwatch() : start_time(clock()) {}
// 	void start() { start_time = clock(); }
// 	long int GetCPUunits() const;
//     void print(std::ostream& outp) const;
// private:
// 	clock_t start_time;
// };
//
// std::ostream& operator<<(std::ostream& outp, const Stopwatch& s);
//
// void estimated_run_time(const mbl::Stopwatch& time, long int N);
//
// class Time_int
// {
// public:
// 	Time_int(long int n) : cnt(n) {}
// 	void operator++(int)
// 	{
// 		cnt++;
// 		if (cnt % 100 == 0 && cnt > 0)
// 			std::cout << "Iteration: " << std::setw(6) << cnt
// 			          << "    " << time << " " << std::endl;
// 	}
// 	bool less(long int N)
// 	{
// 		if (cnt == 10)
// 			mbl::estimated_run_time(time,N/10);
// 		return cnt < N;
// 	}
// 	operator long int(){ return cnt; }
// private:
// 	long int cnt;
// 	Stopwatch time;
// };

// simple Date class, can be used to find weeknr and week of the day.
// see file "ISOwdALG.txt" for more details and "ISO 8601 Week Date"
// for more details
// class Date
// {
// public:
// 	// constructor with year (y), month(m), and day(d)
// 	Date(int y, int m, int d) : yr(y),month(m),day(d) { init(); }
// 	// Get the day of the week: monday = 1, sunday = 7
// 	int GetWeekDay() const { return Weekday; }
// 	// Get name of the day of the week.
// 	std::string GetWeekDayString() const;
// 	// Get Week number
// 	int GetWeekNumber() const { return WeekNumber; }
// 	// Is True if week of the day is saturday or sunday
// 	bool Weekend() const { return Weekday > 5;}
// 	void PrintDate(std::ostream& outp) const;
// 	// print ISO8601 Week date: YearNumber-WeekNumber-Weekday.
// 	void PrintWeekDate(std::ostream& outp) const;
//
// 	friend int compare(const Date& d1, const Date& d2);
// private:
// 	void init();
// 	const static int cumMonth[];
// 	int yr, month, day;
// 	int Weekday, WeekNumber, YearNumber;
// };
//
// int compare(const Date& d1, const Date& d2);
//
// inline bool operator==(const Date& d1, const Date& d2)
// { return (compare(d1,d2) == 0); }
//
// inline bool operator!=(const Date& d1, const Date& d2)
// { return (compare(d1,d2) != 0); }
//
// inline bool operator<(const Date& d1, const Date& d2)
// { return (compare(d1,d2) < 0); }
//
// struct auxPrintWeekDate
// {
// 	auxPrintWeekDate(const Date& d) : date(d) {}
// 	std::ostream& operator()(std::ostream& outp) const
// 		{ date.PrintWeekDate(outp); return outp;}
// 	Date date;
// };
//
// inline std::ostream& operator<<(std::ostream& outp, const Date& date)
// { date.PrintDate(outp); return outp; }
//
// inline auxPrintWeekDate WeekDate(const Date& d)
// { return auxPrintWeekDate(d); }
//
// inline std::ostream& operator<<(std::ostream& outp, const auxPrintWeekDate& aux)
// { return aux(outp); }

// std::istream& skip_header(std::istream &f);
//std::istream& skip_rest_of_line(std::istream& inp);
//std::istream& eatcomment(std::istream& inp);

//void tolower(std::string& a);
//void toupper(std::string& a);

// eenvoudige classe om uitwisseling tussen MFC en rest code goed gescheiden te houden
// class ExchThreadInfo
// {
// public:
// 	ExchThreadInfo() {}
// 	virtual ~ExchThreadInfo() {}
// 	virtual void Update() = 0;
// 	virtual void Message(const char *str) const = 0; // eventueel extra arg: int destination
// };

// bool not_space(char c);
// std::string remove_comment(const std::string& str);

// bool IsNumeric(const std::string& a);
// bool IsNumeric(const std::vector<std::string>& x);
// bool IsNumeric_or_MissingValue(const std::vector<std::string>& x,
// 							   const std::string& miss_val);

// int System(const std::string& command);

// bool FileExists(const std::string& filename);

// std::string make_repetition(int nr);

// std::vector<int> make_range(int begin,int end);

// class MakeLabel
// {
// public:
// 	MakeLabel(const std::string& pre, int width) : pre_(pre),width_(width) {}
// 	std::string operator()(int a);
// private:
// 	std::string pre_;
// 	int width_;
// };
//
}

#endif
