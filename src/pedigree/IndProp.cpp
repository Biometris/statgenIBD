// Martin Boer, Biometris
#include <fstream>
#include <sstream> 
#include <algorithm>

#include "IndProp.h"
#include "misc.h"

using namespace mbl;
using namespace std;


set<string> get_ind_within_fam(const vector<IndProp>& pop)
{
	set<string> x;
	const int N=pop.size();
	for (int i=0;i<N;i++)
	{
		IndProp ind = pop[i];
		if (ind.IsMemberFamily())
			x.insert(ind.GetID());
	}
	return x;
}

set<string> get_parents_families(const vector<IndProp>& pop)
{
	set<string> Pfam;
	const int N=pop.size();
	for (int i=0;i<N;i++)
	{
		IndProp ind = pop[i];
		if (ind.IsMemberFamily())
		{
			Pfam.insert(ind.GetP1());
			Pfam.insert(ind.GetP2());
		}
	}
	return Pfam;
}

set<string> get_ind_coa_file(const vector<IndProp>& pop, 
							 const string& filename_coa)
{
	ifstream inp;
	OpenFile(inp,filename_coa);
	const int N = pop.size();
	set<string> ID_names;
	for (int i=0;i<N;i++)
	{
		IndProp ind = pop[i];
		if (!ind.IsMemberFamily())
		{
			ID_names.insert(ind.GetID());
		}
	}
	
	set<string> set_ind;
	string line;
	while (getline(inp,line))
	{
		if (line.empty()) continue;
		istringstream line_stream(line);
		string ind;
		line_stream >> ind;
		if (ID_names.find(ind) == ID_names.end()) 
			throw mblib_error("ID=" + ind + " not defined in pedigree file");
		set_ind.insert(ind);
	}
	return set_ind;
}

vector<int> get_ndx_set(const vector<IndProp>& pop, const set<string>& setID)
{
	const int N=pop.size();
	vector<int> ndx;
	for (int i=0;i<N;i++)
	{
		if (setID.find(pop[i].GetID()) != setID.end())
			ndx.push_back(i);
	}
	return ndx;
}

bool single_family(const vector<IndProp>& pop)
{
	int npar=0;
	int nfnd=0;
	int nhyb=0;
	int n_3way_and_4way=0;
	for (vector<IndProp>::const_iterator it=pop.begin();it!=pop.end();it++)
	{
		if (it->IsInbredParent())
			npar++;
		else if (it->IsFounder())
			nfnd++;
		else if (it->IsHybrid())
			nhyb++;
		if (it->GetType().find("C") == 0)
			n_3way_and_4way++;
	}
	if (npar != 0)
	{
		if (npar < 2 || npar > 4 || npar-nhyb!=2 || nfnd != 0)
			throw mblib_error("error in check_fam!");
		int nfam=0;
		string fam_name, fam_type, fam_P1, fam_P2;
		for (vector<IndProp>::const_iterator it=pop.begin();it!=pop.end();it++)
		{
			if (it->IsMemberFamily())
			{
				if (nfam == 0)
				{
					fam_name = it->GetFam();
					fam_type = it->GetType();
					fam_P1 = it->GetP1();
					fam_P2 = it->GetP2();
				}
				else
				{
					if (it->GetFam() != fam_name || it->GetType() != fam_type ||
						it->GetP1() != fam_P1 || it->GetP2() != fam_P2)
						throw mblib_error("Error for ID " + it->GetID());
				}
				nfam++;
			}
		}
		return true;
	}
	else
	{
		if (n_3way_and_4way > 0)
			throw mblib_error("Not implemented yet: 3-way or 4-way crosses in multiple populations");
		return false;
	}
}

