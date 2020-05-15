/***********************************************************
*														   *
*   Martin Boer											   *
*   Biometris Wageningen-UR & NWO, grant no. 925-01-001    *
*   July 25, 2006										   *
*														   *
*   July 25, 2006:										   *
*	correction of remove_comment() procedure		       *
*														   *
************************************************************/

#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <limits>
// #include <direct.h>
#include <stdlib.h>
// #include "util.h"
#include "mblexcept.h"

using namespace std;

namespace
{

//void upper_to_lower(char& c) { c = tolower(c); }
//void lower_to_upper(char& c) { c = toupper(c); }

}

// double mbl::mean(const vector<double>& y)
// {
// 	double sum = accumulate(y.begin(),y.end(),0.0);
// 	return sum/y.size();
// }
//
// double mbl::mean(const vector<double>& y, const vector<double>& w)
// {
// 	double sum_w  = 0.0;
// 	double sum_wy = 0.0;
// 	vector<double>::const_iterator iter_y = y.begin();
// 	vector<double>::const_iterator iter_w = w.begin();
// 	for (;iter_y != y.end(); ++iter_y, ++iter_w)
// 	{
// 		sum_wy += (*iter_y)*(*iter_w);
// 		sum_w  += (*iter_w);
// 	}
// 	return sum_wy/sum_w;
// }
//
// double mbl::variance(const vector<double>& y)
// {
// 	double mu = mean(y);
// 	double sum_sqr = 0.0;
// 	vector<double>::const_iterator iter_y = y.begin();
// 	for (;iter_y != y.end(); ++iter_y)
// 		sum_sqr += sqr(*iter_y - mu);
// 	return sum_sqr/y.size();
// }
//
// double mbl::variance(const vector<double>& y, const vector<double>& w)
// {
// 	double mu = mean(y,w);
// 	double sum_sqr = 0.0;
// 	vector<double>::const_iterator iter_y = y.begin();
// 	vector<double>::const_iterator iter_w = w.begin();
// 	for (;iter_y != y.end(); ++iter_y, ++iter_w)
// 		sum_sqr += (*iter_w)*sqr(*iter_y - mu);
// 	return sum_sqr/y.size();
// }
//
// double mbl::covariance(const vector<double>& x, const vector<double>& y)
// {
// 	double mu_x = mean(x);
// 	double mu_y = mean(y);
// 	vector<double>::const_iterator iter_x = x.begin();
// 	vector<double>::const_iterator iter_y = y.begin();
// 	double sum = 0.0;
// 	for (;iter_x != x.end(); ++iter_x, ++iter_y)
// 		sum += (*iter_x - mu_x)*(*iter_y - mu_y);
// 	return sum/x.size();
// }

// double mbl::sum_of_sqr(const vector<double>& x)
// {
// 	double sum_sqr = 0.0;
// 	vector<double>::const_iterator iter_x = x.begin();
// 	for (;iter_x != x.end(); ++iter_x)
// 		sum_sqr += sqr(*iter_x);
// 	return sum_sqr;
// }
//
// double mbl::sum_of_sqr(const vector<double>& x, const vector<double>& w)
// {
// 	double sum_sqr = 0.0;
// 	vector<double>::const_iterator iter_x = x.begin();
// 	vector<double>::const_iterator iter_w = w.begin();
// 	for (;iter_x != x.end(); ++iter_x, ++iter_w)
// 		sum_sqr += (*iter_w)*sqr(*iter_x);
// 	return sum_sqr;
// }

// void mbl::Stopwatch::print(ostream& outp) const
// {
// 	clock_t end = clock();
// 	long int cpu_units= (long int)(end - start_time);
// 	long int cpu= cpu_units/(long int)CLOCKS_PER_SEC;
// 	long int sec= cpu%60;
// 	long int min= cpu/60;
// 	long int hours= min/60;
// 	if (hours)  min%= 60;
// 	char old_fill = outp.fill();
// 	outp << setw(2) << setfill('0') << hours << ":"
// 		 << setw(2) << setfill('0') << min	 << ":"
// 		 << setw(2) << setfill('0') << sec
// 		 << "  ("   << cpu_units << " units)";
// 	outp << setfill(old_fill);
// }
//
// long int mbl::Stopwatch::GetCPUunits() const
// {
// 	clock_t end = clock();
// 	long int cpu_units= (long int)(end - start_time);
// 	return cpu_units;
// }
//
// ostream& mbl::operator<<(ostream& outp, const Stopwatch& s)
// {
// 	s.print(outp);
// 	return outp;
// }

// const int mbl::Date::cumMonth[] = {0,31,59,90,120,151,181,212,243,273,304,334};
//
// bool leapyear(int yr) { return (yr%4==0 && yr % 100 !=0) || (yr % 400 ==0);}
//
// //mbl::Date::Date(int y, int m, int d) : yr(y),month(m),day(d)
//
// void mbl::Date::init()
// {
// 	int daynr_year = day + cumMonth[month-1];
// 	if (leapyear(yr) && month > 1)
// 		daynr_year += 1;
// 	// find the Jan1Weekday for yr (monday = 1, sunday = 7)
// 	int YY = (yr-1) % 100;
// 	int C = (yr-1) - YY;
// 	int G = YY + YY/4;
// 	int Jan1Weekday = 1 + (((((C/100) % 4) * 5) + G) % 7);
// 	int H = daynr_year + (Jan1Weekday-1);
// 	Weekday = 1+(H-1)%7;
// 	int YearNumber,WeekNumber;
// 	if (daynr_year <= 8-Jan1Weekday && Jan1Weekday > 4)
// 	{
// 		YearNumber = yr-1;
// 		if (Jan1Weekday==5 || (Jan1Weekday=6 && leapyear(yr-1)))
// 			WeekNumber = 53;
// 		else
// 			WeekNumber = 52;
// 	}
// 	else
// 		YearNumber = yr;
// 	if (YearNumber == yr)
// 	{
// 		int I = leapyear(yr) ? 366 : 365;
// 		if (I - daynr_year < 4-Weekday)
// 		{
// 			YearNumber = yr+1;
// 			WeekNumber = 1;
// 		}
// 	}
// 	if (YearNumber == yr)
// 	{
// 		int J = daynr_year + (7-Weekday) + (Jan1Weekday - 1);
// 		WeekNumber = J/7;
// 		if (Jan1Weekday > 4)
// 			WeekNumber -= 1;
// 	}
// }
//
// void mbl::Date::PrintDate(std::ostream& outp) const
// {
// 	outp << setfill('0')
// 		 << setw(2) << day << "-"
// 		 << setw(2) << month << "-"
// 		 << yr << setfill(' ');
// }
//
// void mbl::Date::PrintWeekDate(ostream& outp) const
// {
// 	outp << YearNumber << "-W"
// 		 << setw(2) << setfill('0') <<  WeekNumber << "-" << Weekday << setfill(' ');
// }
//
// string mbl::Date::GetWeekDayString() const
// {
//     switch (Weekday)
//     {
//       case 1: return "mon";
// 	  case 2: return "tue";
// 	  case 3: return "wed";
// 	  case 4: return "thu";
// 	  case 5: return "fri";
// 	  case 6: return "sat";
// 	  case 7: return "sun";
// 	}
// 	return "";
// }
//
// int mbl::compare(const Date& d1, const Date& d2)
// {
// 	if (d1.yr < d2.yr) return -1;
// 	if (d1.yr > d2.yr) return  1;
//
// 	if (d1.month < d2.month) return -1;
// 	if (d1.month > d2.month) return  1;
//
// 	if (d1.day < d2.day) return -1;
// 	if (d1.day > d2.day) return  1;
//
// 	return 0;
// }

// istream& mbl::skip_header(istream &inp)
// {
// 	while (inp.peek() == '#')
// 		skip_rest_of_line(inp);
// 	return inp;
// }

// bool mbl::not_space(char c)
// {
// 	return !::isspace((unsigned) c);
// }
//
// string mbl::remove_comment(const string& str)
// {
// 	string::const_iterator i,j,k;
// 	i = str.begin();
// 	i = find_if(i,str.end(),not_space);
//
// 	if (i == str.end() || (*i == '#'))
// 		return string();
//
// 	j = find(i,str.end(),';');
//
// 	// added 25 july 2006
// 	if (j == str.begin())
// 		return string();
// 	// end added 25 july 2006
//
// 	for (k=j-1;k>=i;k--)
// 		if (not_space(*k))
// 			break;
//
//
// 	return string(i,k+1);
// }

// bool mbl::IsNumeric(const string& a)
// {
// 	istringstream inp(a);
// 	double val;
// 	inp >> val;
// 	eatcomment(inp);
// 	return (!inp.bad() && inp.eof());
// }
//
// bool mbl::IsNumeric(const vector<string>& x)
// {
// 	for (vector<string>::const_iterator it = x.begin(); it != x.end(); ++it)
// 		if (IsNumeric(*it) == false)
// 			return false;
// 	return true;
// }
//
// bool mbl::IsNumeric_or_MissingValue(const vector<string>& x, const string& miss_val)
// {
// 	bool all_missing_values = true;
// 	for (vector<string>::const_iterator it = x.begin(); it != x.end(); ++it)
// 	{
// 		if (*it != miss_val)
// 		{
// 			if (IsNumeric(*it))
// 				all_missing_values = false;
// 			else
// 				return false;
// 		}
// 	}
// 	return !all_missing_values;
// }

// int mbl::System(const string& command) { return system(command.c_str()); }


// bool mbl::FileExists(const string& filename)
// {
// 	ifstream inp(filename.c_str());
// 	return (inp) ? true : false;
// }

// string mbl::make_repetition(int nr)
// {
// 	string nr_str = itostr(nr+1);
// 	if (nr+1 < 10)    return string("rep000" + nr_str);
// 	if (nr+1 < 100)   return string("rep00"  + nr_str);
// 	if (nr+1 < 1000)  return string("rep0"   + nr_str);
// 	if (nr+1 < 10000) return string("rep"    + nr_str);
//
// 	mblib_error("number of repetitions >= 10000");
// 	return ""; // dummy
// }

// void mbl::estimated_run_time(const Stopwatch& time, long int N)
// {
// 	long int cpu_units = time.GetCPUunits()*N;
// 	long int cpu= cpu_units/(long int)CLOCKS_PER_SEC;
// 	long int sec= cpu%60;
// 	long int min= cpu/60;
// 	long int hours= min/60;
// 	if (hours)  min%= 60;
// 	char old_fill = cout.fill();
// 	cout << "Estimated time needed for calculations: ";
// 	cout << setw(2) << setfill('0') << hours << ":"
// 		 << setw(2) << setfill('0') << min	 << ":"
// 		 << setw(2) << setfill('0') << sec  << endl;
// 	cout << setfill(old_fill);
// }

// [begin,end)
// vector<int> mbl::make_range(int begin,int end)
// {
// 	vector<int> x;
// 	for (int i=begin;i!=end;i++)
// 		x.push_back(i);
// 	return x;
// }

// string mbl::MakeLabel::operator()(int a)
// {
// 	string str_a = stringify(a+1);
// 	if ((int)str_a.length() > width_)
// 		throw mbl::mblib_error("MakeLabel");
// 	ostringstream o;
// 	o << pre_ << std::setfill('0') << setw(width_) << str_a;
// 	return o.str();
// }


/*
string mbl::GetEnv(const string& var_name)
{
	char *path = getenv(var_name.c_str());
	if (!path)
		throw mblib_error("Cannot find variable " + var_name + " in environment table");
	return string(path);
}
*/


