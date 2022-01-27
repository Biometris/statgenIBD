#include <algorithm>
#include <numeric>
#include <fstream>
#include <iostream> // for setfill
#include <iomanip>
#include <sstream>
#include <cfloat>
#include <float.h> // for DBL_MIN
#include <limits> // for numeric_limits

#include "misc.h"
#include "ibdexcept.h"
#include "convert.h"

using namespace ibd;
using namespace std;

double ibd::round(double x,
                  int precision)
{
	double fac = pow10(precision);
	return floor(fac*x+0.5)/fac;
}

double ibd::mean(const vector<double>& y)
{
	double sum = accumulate(y.begin(),y.end(),0.0);
	return sum/y.size();
}

double ibd::mean(const vector<double>& y, const vector<double>& w)
{
	double sum_w  = 0.0;
	double sum_wy = 0.0;
	vector<double>::const_iterator iter_y = y.begin();
	vector<double>::const_iterator iter_w = w.begin();
	for (;iter_y != y.end(); ++iter_y, ++iter_w)
	{
		sum_wy += (*iter_y)*(*iter_w);
		sum_w  += (*iter_w);
	}
	return sum_wy/sum_w;
}

double ibd::variance(const vector<double>& y)
{
	double mu = mean(y);
	double sum_sqr = 0.0;
	vector<double>::const_iterator iter_y = y.begin();
	for (;iter_y != y.end(); ++iter_y)
		sum_sqr += sqr(*iter_y - mu);
	return sum_sqr/y.size();
}

double ibd::variance(const vector<double>& y, const vector<double>& w)
{
	double mu = mean(y,w);
	double sum_sqr = 0.0;
	vector<double>::const_iterator iter_y = y.begin();
	vector<double>::const_iterator iter_w = w.begin();
	for (;iter_y != y.end(); ++iter_y, ++iter_w)
		sum_sqr += (*iter_w)*sqr(*iter_y - mu);
	return sum_sqr/y.size();
}

// conversion from int to string
string ibd::itostr(int a)
{
	return stringify(a);
}

void upper_to_lower(char& c) { c = tolower(c); }
void lower_to_upper(char& c) { c = toupper(c); }

void ibd::tolower(string& a)
{
	for_each(a.begin(),a.end(),upper_to_lower);
}

void ibd::toupper(string& a)
{
	for_each(a.begin(),a.end(),lower_to_upper);
}

void ibd::make_conditional(vector<double>& p)
{
	double sum = accumulate(p.begin(),p.end(),0.0);
	if (sum <= DBL_MIN)
		throw ibd_error("function make_conditional()");
	for (vector<double>::iterator it = p.begin(); it != p.end(); ++it)
		*it/=sum;
}

vector<double> ibd::elem_prod(const vector<double>& a, const vector<double>& b)
{
	size_t dim = a.size();
	vector<double> result(dim);
	for (unsigned int i=0;i<dim;i++)
		result[i] = a[i]*b[i];
	return result;
}

istream& ibd::skip_header(istream &inp)
{
	while (inp.peek() == '#')
		skip_rest_of_line(inp);
	return inp;
}

istream& ibd::skip_rest_of_line(istream& inp)
{
	inp.ignore(numeric_limits<int>::max(),'\n');
	return inp;
}

istream& ibd::skip_lines(istream& inp, int n)
{
	for (int i=0;i<n;i++)
		skip_rest_of_line(inp);
	return inp;
}

istream& ibd::eatcomment(istream& inp)
{
    char c;
    while (inp.get(c))
	{
		if (c == ';' || c == '#')
			skip_rest_of_line(inp);
		else if (!isspace(c))
		{
			inp.putback(c);
			return inp;
		}
	}
	return inp;
}

bool ibd::not_space(char c)
{
	return !::isspace((unsigned) c);
}

string ibd::remove_comment(const string& str)
{
	string::const_iterator i,j,k;
	i = str.begin();
	i = find_if(i,str.end(),not_space);

	if (i == str.end() || (*i == '#'))
		return string();

	j = find(i,str.end(),';');

	// added 25 july 2006
	if (j == str.begin())
		return string();
	// end added 25 july 2006

	for (k=j-1;k>=i;k--)
		if (not_space(*k))
			break;


	return string(i,k+1);
}

void ibd::OpenFile(ofstream& outp, string filename)
{
	outp.open(filename.c_str());
	if (!outp)
		throw ibd_error("Cannot open file " + filename);
	outp.setf(ios::fixed, ios::floatfield);
}

void ibd::OpenFile(ifstream& inp, string filename)
{
	inp.open(filename.c_str());
	if (!inp)
		throw ibd_error("Cannot open file " + filename);
}

unsigned int ibd::pow2(int n)
{
	unsigned int x=1;
	return (x << n);
}

string ibd::MakeLabel::operator()(int a)
{
	string str_a = stringify(a+1);
	if ((int)str_a.length() > width_)
		throw ibd::ibd_error("MakeLabel");
	ostringstream o;
	o << pre_ << std::setfill('0') << setw(width_) << str_a;
	return o.str();
}
