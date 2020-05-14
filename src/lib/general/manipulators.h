// Martin Boer, 18 october 2004
// istream & ostream manipulators
#ifndef MANIPULATORS_HEADER
#define MANIPULATORS_HEADER

//#pragma warning(disable:4786)   // disable C4786 warning

#include "util.h"

namespace mbl
{

// istream function object (for manipulators with 1 argument of type T)
template <class T>
class istream_fo1
{
public:
	istream_fo1(std::istream& (*f)(std::istream& , T), int v) : fun(f), val(v) {}
	std::istream& operator()(std::istream& i) const { return (*fun)(i,val);}
private:
	std::istream& (*fun)(std::istream&, T);
	T val;
};

template<class T>
std::istream& operator>>(std::istream& inp, istream_fo1<T>& fo) { return fo(inp); }

std::istream& skip_fields(std::istream& inp, int n);
std::istream& skip_lines(std::istream& inp, int n);

//istream_fo1<int> skipf(int n = 1);
//istream_fo1<int> skipl(int n = 1);
inline istream_fo1<int> skipf(int n = 1) { return istream_fo1<int>(skip_fields, n);}
inline istream_fo1<int> skipl(int n = 1) { return istream_fo1<int>(skip_lines, n);}

}

#endif

