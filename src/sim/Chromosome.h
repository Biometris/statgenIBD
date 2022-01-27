#ifndef CHROMOSOME_HEADER
#define CHROMOSOME_HEADER

#include <vector>
#include <iostream>
#include "Genotype.h"
#include "util_genetics.h"
#include "binaryfile.h"

namespace ibd
{

class ChromosomePair;
class Genome;

struct Segment
{
	Segment() {} // 12 oktober 2000: Borland c++ needs this constr.
	Segment(alleletype a, double r) : allele(a), right(r) {}
	alleletype allele;
	double right; // right end of segment
};

inline oBinFile& operator<<(oBinFile& os, const Segment& x)
{ os.write((char *) &x,sizeof(Segment)); return os; }

inline iBinFile& operator>>(iBinFile& is, Segment& x)
{ is.read((char *) &x,sizeof(Segment)); return is;}


class Chromosome  
{
public:
	typedef std::vector<Segment>::const_iterator const_iterator;
	typedef std::vector<Segment>::iterator iterator;

	Chromosome() {}
	Chromosome(double length, alleletype allele) : segments(1,Segment(allele,length)) {}
	alleletype GetGenotype(double cM) const;
	int GetNrSegments() const { return segments.size(); }
	double Length() const { return segments.back().right; } 
	void print(std::ostream&) const ;
	std::vector<Segment> GetSegments() const { return segments; }
private:
	void add_segment(const Segment& s) { segments.push_back(s); }

	std::vector<Segment> segments;

	friend bool eq(const Chromosome&, const Chromosome&);
	friend ChromosomePair;
	friend Genome;

	friend oBinFile& operator<<(oBinFile& os, const Chromosome& x);
	friend iBinFile& operator>>(iBinFile& is, Chromosome& x);

};

bool eq(const Chromosome&, const Chromosome&);

inline bool operator==(const Chromosome& chr1, const Chromosome& chr2)
{ return eq(chr1,chr2); }

inline bool operator!=(const Chromosome& chr1, const Chromosome& chr2)
{ return !eq(chr1,chr2); }

oBinFile& operator<<(oBinFile& os, const Chromosome& x);
iBinFile& operator>>(iBinFile& is, Chromosome& x);

}

#endif 
