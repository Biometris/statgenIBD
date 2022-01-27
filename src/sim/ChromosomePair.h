#ifndef CHROMOSOMEPAIR_HEADER
#define CHROMOSOMEPAIR_HEADER

#include "Chromosome.h"
#include "Genotype.h"
#include "util_genetics.h"
#include "binaryfile.h"

namespace ibd
{

class Genome;

class ChromosomePair  
{
public:
	ChromosomePair() {}
	ChromosomePair(double l, alleletype a, alleletype b) : chr_a(l,a), chr_b(l,b) {}
	Chromosome Meiosis() const;
	void print(std::ostream&) const;
	Genotype GetGenotype(double cM) const 
		{ return Genotype(chr_a.GetGenotype(cM),chr_b.GetGenotype(cM)); }
	int GetNrSegments () const 
		{ return chr_a.GetNrSegments() + chr_b.GetNrSegments(); } 
	bool IsHomozygote() const { return (chr_a == chr_b); }
private:
	Chromosome chr_a;
	Chromosome chr_b;

	friend bool eq(const ChromosomePair& a, const ChromosomePair& b);
	friend Genome Cross(const Genome& a, const Genome& b);
	friend Genome DoubledHaploid(const Genome& a);
	friend Genome;
	friend oBinFile& operator<<(oBinFile& os, const ChromosomePair& x);
	friend iBinFile& operator>>(iBinFile& is, ChromosomePair& x);
};

bool eq(const ChromosomePair& a, const ChromosomePair& b);

inline bool operator==(const ChromosomePair& a, const ChromosomePair& b) 
{ return eq(a,b); }

inline bool operator!=(const ChromosomePair& a, const ChromosomePair& b) 
{ return !eq(a,b); }

oBinFile& operator<<(oBinFile& os, const ChromosomePair& x); 
iBinFile& operator>>(iBinFile& is, ChromosomePair& x);       


}

#endif
