// Martin Boer, 18 october 2004
// istream & ostream manipulators
#ifndef MANIPULATORS_HEADER
#define MANIPULATORS_HEADER

#include <iostream>

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

std::istream& skip_lines(std::istream& inp, int n);

inline istream_fo1<int> skipl(int n = 1) { return istream_fo1<int>(skip_lines, n);}

}

#endif

