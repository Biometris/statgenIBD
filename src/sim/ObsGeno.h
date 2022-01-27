#ifndef OBSGENO_HEADER
#define OBSGENO_HEADER

#include <vector>
#include <iostream>
#include "util_genetics.h"
#include "Genotype.h"
#include "SmartPtr.h"

namespace ibd
{

class ObsGeno 
{
public:
	typedef std::vector<Genotype>::const_iterator const_iterator;
	typedef std::vector<Genotype>::size_type size_type;

	ObsGeno() : h_data(new std::vector<Genotype>) {}
	ObsGeno(const Genotype& g1);
	ObsGeno(const Genotype& g1, const Genotype& g2);
	ObsGeno(const Genotype& g1, const Genotype& g2, const Genotype& g3);
	ObsGeno(const Genotype& g1, const Genotype& g2, const Genotype& g3, const Genotype& g4);

	const_iterator begin() const { return h_data->begin(); }
	const_iterator end() const { return h_data->end(); }
	size_type size() const { return h_data->size(); }
	bool empty() const { return h_data->empty(); }
	const Genotype& operator[](size_t i) const {return (*h_data)[i]; }
	void print(std::ostream& outp) const;
	friend bool eq2(const ObsGeno& a, const ObsGeno& b);
	friend ObsGeno intersection(const ObsGeno& g1, const ObsGeno& g2);

private:
	ObsGeno(const std::vector<Genotype>& g);
	SmartPtr< std::vector<Genotype> > h_data; // pointer to data

};

inline bool eq2(const ObsGeno& a, const ObsGeno& b)
{ return (*a.h_data == *b.h_data);}

inline bool operator==(const ObsGeno& a, const ObsGeno& b) 
{ return eq2(a,b); }

inline bool operator!=(const ObsGeno& a, const ObsGeno& b)
{ return !eq2(a,b); }

std::ostream& operator<<(std::ostream& outp, const ObsGeno& g);

ObsGeno intersection(const ObsGeno& g1, const ObsGeno& g2);

}

#endif

