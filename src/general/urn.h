#ifndef URN_HEADER
#define URN_HEADER

#include <vector>
#include "ibdexcept.h"
#include "Random.h"

namespace ibd
{

template <class T>
class Urn
{
public:
	Urn() : n(0) {}                                      // empty urn
	Urn(const std::vector<T>& x) : y(x), n(x.size()) { } // constructor
	T random_draw();									 // draw a random ball
 	int total_nr_balls() const { return y.size(); }      // total nr of balls in experiment
	int nr_balls_in_urn() const { return n; }            // number of balls in urn
	void switch_balls();                                 // switch balls in and out of urn
	void draw(const T& i);						         // draw element i
	void draw(const T& i, int ndx);                      // see comment in urn.cpp
	void replace(const T& i);							 // replace element i
	void replace_all() { n = y.size(); }			     // replace all the balls 
private:
	std::vector<T> y;								     // vector with all the "balls"
	int n;												 // number of "balls" in "urn"
};

template <class T>
T Urn<T>::random_draw()
{
	if (n==0) 
		throw ibd_error("Urn::random_draw() : empty urn");
	int r = randint(n);
	std::swap(y[r],y[--n]);
	return y[n];
}

/*

template <class T>
void Urn<T>::draw(const T& i)
{
	typedef std::vector<T>::iterator vec_iter;
	vec_iter begin = y.begin();
	vec_iter end   = begin + n;
	vec_iter it    = find(begin, end, i);
	if (it == end) 
		throw ibd_error("Urn::draw(const T& i)");
	std::swap(*it,y[--n]);
}

// draw element i, using the initial guess "pos" of the position of the element
// for example, if y has the following values:
// A,B,C,D,E,F,G 
// and n=5 (number of balls in urn)
// then draw(A,2) will search in the following order: 2,3,4,0,1 
template <class T>
void Urn<T>::draw(const T& i, int pos)
{
	if (pos >= n || pos < 0) 
		throw ibd_error("Urn::draw(const T& i, int ndx)");
	typedef std::vector<T>::iterator vec_iter;
	vec_iter begin           = y.begin();
	vec_iter pos_first_guess = begin + pos;
	vec_iter end             = begin + n;

	vec_iter it  = find(pos_first_guess, end, i);
	if (it == end)						
	{
		it = find(begin,pos_first_guess,i); 
		if (it == pos_first_guess)
			throw ibd_error("Urn::draw(const T& i, int ndx)");
	}
	std::swap(*it,y[--n]);
}

template <class T>
void Urn<T>::replace(const T& i)
{
	vector<T>::iterator it = find (y.begin() + n, y.end(), i); 
	if (it == y.end()) 
		throw ibd_error("Urn::replace(const T& i)");
	std::swap(*it,y[n++]); 
}

template <class T>
void Urn<T>::switch_balls()
{
    reverse(y.begin(), y.end());
	n = y.size() - n;
}

int test1_Urn();
int test2_Urn();

*/

}

#endif