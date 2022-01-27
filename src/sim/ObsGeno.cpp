#include "ObsGeno.h"
#include <algorithm>
	
using namespace std;
using namespace ibd;

ObsGeno::ObsGeno(const Genotype& g1) 
		: h_data(new vector<Genotype>(1))
{ 
	(*h_data)[0] = g1; 
}

ObsGeno::ObsGeno(const Genotype& g1, const Genotype& g2) 
		: h_data(new vector<Genotype>(2))
{ 
	(*h_data)[0] = g1; 
	(*h_data)[1] = g2; 
	sort(h_data->begin(), h_data->end());
}
	
ObsGeno::ObsGeno(const Genotype& g1, const Genotype& g2, const Genotype& g3)
		: h_data(new vector<Genotype>(3))
{ 
	(*h_data)[0] = g1; 
	(*h_data)[1] = g2; 
	(*h_data)[2] = g3; 
	sort(h_data->begin(), h_data->end());
}

ObsGeno::ObsGeno(const Genotype& g1, const Genotype& g2, 
				 const Genotype& g3, const Genotype& g4)
		: h_data(new vector<Genotype>(4))
{
	(*h_data)[0] = g1; 
	(*h_data)[1] = g2; 
	(*h_data)[2] = g3; 
	(*h_data)[3] = g4;
	sort(h_data->begin(), h_data->end());
}	

ObsGeno::ObsGeno(const vector<Genotype>& g) 
		: h_data(new vector<Genotype>(g))
{ 
	sort(h_data->begin(), h_data->end());
}



void ObsGeno::print(ostream& outp) const
{
	const int nrgeno = h_data->size();
	outp << "{";
	vector<Genotype>::const_iterator iter = h_data->begin();
	for (int i=0;i<nrgeno;i++)
	{
		outp << *iter++;
		if (i < nrgeno-1)
			outp << ",";
	}
	outp << "}";
}

ostream& ibd::operator<<(ostream& outp, const ObsGeno& g)
{
	int str_length = g.empty() ? 2 : 3*g.size() + 1;
	int nr_whitespaces = outp.width() - str_length;
	outp.width(0);
	for (int i=0;i<nr_whitespaces;i++)
		outp << ' ';
	g.print(outp);
	return outp;
}

ObsGeno ibd::intersection(const ObsGeno& g1, const ObsGeno& g2)
{
	const vector<Genotype>& v1 = *g1.h_data;
	const vector<Genotype>& v2 = *g2.h_data;
	vector<Genotype> v3;

	for (vector<Genotype>::const_iterator it = v1.begin(); it != v1.end(); ++it)
		if (find(v2.begin(),v2.end(),*it) != v2.end())
			v3.push_back(*it);

	return ObsGeno(v3);
}


