// Martin Boer, Biometris
#include <map>
#include <fstream>
#include <sstream> 
#include "misc.h"
#include "pedigree.h"
#include "mblexcept.h"

using namespace mbl;
using namespace std;

const int& ParentsInd::operator[](int x) const
{
	if (x==0) return this->first;
	if (x==1) return this->second;
	throw mblib_error("error in Parents::operator[]");
	return this->first; // dummy
}

Pedigree::Pedigree(const std::vector<IndProp>& pop) 
{
	Nfnd = 0;
	const int N = pop.size();
	map<string,int> ndx;
	for (int i=0;i<N;i++)
		ndx[pop[i].GetID()] = i;

	for (int i=0;i<N;i++)
	{
		const IndProp& ind = pop[i];
		ParentsInd parents;
		if (ind.IsFounder())
			Nfnd++;
		else
		{
			int ndx_P1 = ndx[ind.GetP1()];
			int ndx_P2 = ndx[ind.GetP2()];
			parents = ParentsInd(ndx_P1,ndx_P2);
		}
		this->push_back(parents);
	}
}

Pedigree::Pedigree(const std::vector<ParentsInd>& vecparents)
	: vector<ParentsInd>(vecparents), Nfnd(0)
{
	const int N = vecparents.size();
	for (int i=0;i<N;i++)
		if (vecparents[i].IsFounder())
			Nfnd++;
}

std::vector<int> Pedigree::GetNdxFnd() const
{
	int N = this->size();
	vector<int> x;
	for (int i=0;i<N;i++)
	{
		const ParentsInd& parents = (*this)[i];
		if (parents.IsFounder())
			x.push_back(i);
	}
	return x;
}

