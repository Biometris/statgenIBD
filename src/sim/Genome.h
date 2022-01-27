#ifndef GENOME_HEADER
#define GENOME_HEADER

#include <vector>
#include <map>
#include "ChromosomePair.h"
#include "Loc.h"
#include "Genotype.h"
#include "binaryfile.h"

namespace ibd 
{

enum Parent {Father, Mother};

class Genome  
{
public:
	typedef std::vector<ChromosomePair>::const_iterator const_iterator;
	typedef std::vector<ChromosomePair>::iterator iterator;
	typedef std::vector<ChromosomePair>::size_type size_type;

	Genome() {}
	Genome(const std::vector<double>& , alleletype , alleletype );
	Genome(int, double, alleletype, alleletype);
	Genome(const std::vector<double>& , const Genotype& g);
	Genome(int, double, const Genotype& g);
	Genotype GetGenotype(int chr_nr, double cM) const;
	Genotype GetGenotype(const Locus& loc) const;
	int GetNrChromosomePairs() const { return chr_pair.size(); }
	std::vector<Genotype> GetGenotype(const LinkageMap&) const;
	bool IsHomozygote() const;
	int GetNrSegments() const;
	void print(std::ostream&) const;

	std::vector<alleletype> GetGamete(const LinkageMap&, Parent parent) const; // added 25 june 2004
	std::multimap<alleletype,double> GetLengthFragments(int chr, Parent parent) const;  // added 27 june 2004
	std::vector<Segment> GetSegments(int chr, Parent parent) const; // added 27 june 2004
private:
	std::vector<ChromosomePair> chr_pair;

	friend Genome Cross(const Genome& a, const Genome& b);
	friend Genome DoubledHaploid(const Genome& a);
	friend bool eq(const Genome& a, const Genome& b);

	friend oBinFile& operator<<(oBinFile& os, const Genome& x);
	friend iBinFile& operator>>(iBinFile& is, Genome& x);

};

Genome Cross(const Genome& a, const Genome& b);
Genome DoubledHaploid(const Genome& a);
Genome Selfing(const Genome& x, int Ngen);

bool eq(const Genome& a, const Genome& b);

inline bool operator==(const Genome& a, const Genome& b) { return eq(a,b); }
inline bool operator!=(const Genome& a, const Genome& b) { return !eq(a,b); }

inline Genome operator*(const Genome& a, const Genome& b) { return Cross(a,b); }

std::ostream& operator<<(std::ostream& outp, const Genome& g);

oBinFile& operator<<(oBinFile& os, const Genome& x);
iBinFile& operator>>(iBinFile& is, Genome& x); 

}

#endif 
