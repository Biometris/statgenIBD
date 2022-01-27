#ifndef GENOTYPE_HEADER
#define GENOTYPE_HEADER

#include <vector>
#include <iostream>

namespace ibd 
{

typedef unsigned char alleletype;

class Genotype
{
public:
	Genotype(){}
	Genotype(alleletype a, alleletype b) : x(a), y(b) { if (x > y) std::swap(x,y);}
	Genotype(const Genotype& g) : x(g.x), y(g.y) {}
	alleletype First() const { return x; }
	alleletype Second() const { return y; }
	void print(std::ostream& f) const { f << x << y; }
private:
	alleletype x;
	alleletype y;
    friend int compare(const Genotype& g1, const Genotype& g2);
};

int compare(const Genotype& g1, const Genotype& g2);

inline bool operator==(const Genotype& g1, const Genotype& g2)
{ return (compare(g1,g2) == 0); }

inline bool operator!=(const Genotype& g1, const Genotype& g2)
{ return (compare(g1,g2) != 0); }

inline bool operator<(const Genotype& g1, const Genotype& g2)
{ return (compare(g1,g2) < 0); }

std::ostream& operator<<(std::ostream&, const Genotype&);

}

#endif
