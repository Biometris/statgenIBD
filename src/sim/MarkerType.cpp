#include "MarkerType.h"
#include "urn.h"
#include "convert.h"

using namespace ibd;
using namespace std;

MarkerType::MarkerType(int Nfnd, int Nalleles, double perc_missing)
	: percentage_missing(perc_missing)
{
	if (Nalleles > 0)
	{
		for (int i=0;i<Nfnd;i++)
		{
			int r = randint(1,Nalleles);
			table.push_back(stringify(r)); // randint(1,Nalleles);
		}
	}
	else // if Nalleles <= 0: Generate complete information
	{
		vector<string> x;
		for (int i=0;i<Nfnd;i++)
			x.push_back(stringify(i+1));
		Urn<string> urn(x);
		for (int i=0;i<Nfnd;i++)
			table.push_back(urn.random_draw());
	}
}

string MarkerType::operator()(Genotype g) const 
{	
	double x = randuniform();
	if (x < percentage_missing) return "*";
	int a = g.First();
	int b = g.Second();
	if (GLOBAL_FLAPJACK)
	{
		string c1 = table[a];
		string c2 = table[b];
		if (c1 == c2) 
			return c1;
		if (c1 > c2) 
			swap(c1,c2);
		string score = c1 + "/" + c2;
		return score;
	}
	else
	{
		string score = "(" + table[a] + "," + table[b] + ")";
		return score;
	}
}

// June 7, 2009: For the moment we assume that there are no missing data:
string MarkerType::First(ibd::Genotype g) const
{
	return table[g.First()];
}

// June 7, 2009: For the moment we assume that there are no missing data:
string MarkerType::Second(ibd::Genotype g) const
{
	return table[g.Second()];
}


vector<MarkerType> generate_markertypes(int nloc,int nfnd, int nalleles, 
										  double percentage_missing)
{
	vector<MarkerType> result;
	for (int i=0;i<nloc;i++)
		result.push_back(MarkerType(nfnd,nalleles,percentage_missing));
	return result;
}

