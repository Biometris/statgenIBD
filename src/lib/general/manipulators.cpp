#include "manipulators.h"

using namespace std;

istream& mbl::skip_fields(istream& inp, int n)
{
	string tmp;
	for (int i=0;i<n;i++)
		inp >> tmp;
	return inp;
}

istream& mbl::skip_lines(istream& inp, int n)
{
	for (int i=0;i<n;i++)
		skip_rest_of_line(inp);
	return inp;
}
