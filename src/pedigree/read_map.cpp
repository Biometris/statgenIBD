// Martin Boer, Biometris
#include <list>
#include <set>
#include <map>
#include <numeric>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>

// library file
#include "Args.h"

#include "read_map.h"

using namespace mbl;
using namespace std;

LinkageMap select_chr(const LinkageMap& markermap, int sel_chr)
{
	const int M = markermap.size();
	LinkageMap sel_markermap;
	for (int m=0;m<M;m++)
	{
		Locus loc = markermap[m];
		if (loc.GetChr() == sel_chr)
			sel_markermap.push_back(loc);
	}
	return sel_markermap;
}

LinkageMap adjust_markermap(const LinkageMap& markermap)
{
	const int M = markermap.size();
	const double eps = 1.0e-6;
	const double delta = 0.001;
	LinkageMap new_markermap;
	for (int m=0;m<M;m++)
	{
		Locus loc = markermap[m];
		if (!new_markermap.empty())
		{
			if (loc.GetChr() == new_markermap.back().GetChr())
			{
				double dist = loc.GetPosition() - new_markermap.back().GetPosition();
				if (dist < eps)
				{
					string n = loc.GetName();
					double p = new_markermap.back().GetPosition();
					int c = loc.GetChr();
					loc = Locus(c,p+delta,n);
				}
			}
		}
		new_markermap.push_back(loc);
	}
	return new_markermap;
}

LinkageMap reduce_markermap(const LinkageMap& markermap, const vector<string>& markers)
{
	const int M = markermap.size();
	LinkageMap new_markermap;
	for (int m=0;m<M;m++)
	{
		Locus loc = markermap[m];
		if (((find(markers.begin(),markers.end(),loc.GetName()) != markers.end())
			|| (loc.GetName() == EVAL_POS) || loc.GetName() == EXTR_POS))
		{
			new_markermap.push_back(loc);
		}
	}
	return new_markermap;
}

void print_marker_warnings(const LinkageMap& markermap, const vector<string>& markers)
{
	const int M = markermap.size();
	vector<string> ignored_markers;
	for (vector<string>::const_iterator it = markers.begin();it!=markers.end();it++)
	{
		string name = *it;
		bool marker_found = false;
		for (int m=0;m<M;m++)
		{
			if (markermap[m].GetName() == name)
			{
				marker_found = true;
				break;
			}
		}
		if (!marker_found)
			ignored_markers.push_back(name);
	}
	if (!ignored_markers.empty())
	{
		cout << endl << "WARNING: " << endl
			 << "The following markers will be ignored (not defined in the map file): " << endl;
		for (vector<string>::const_iterator it=ignored_markers.begin();it!=ignored_markers.end();it++)
			cout << *it << endl;
		cout << endl << endl;
	}

}

LinkageMap read_map_file(const string& filename)
{
	LinkageMap markermap;
	ifstream inp(filename.c_str());
	if (!inp)
		throw mblib_error("Cannot read file" + filename);
	string str;
	//int chr_nr = -1;
	while (inp)
	{
	  double pos;
	  string name;
	  int chr;
	  inp >> name >> chr >> pos;
	  Locus loc(chr,pos,name);
	  markermap.push_back(loc);

		//double pos;
		//inp >> eatcomment >> str;
		//tolower(str);
		//if (str == "group" || str == "chrom") // new chromosome
		//{
		//	inp >> chr_nr;
		//}
		//else
		//{
		//	if (chr_nr <= 0)
		//		throw mblib_error("chromosome number not defined");
		//	inp >> pos;
		//	Locus loc(chr_nr,pos,str);
		//	markermap.push_back(loc);
		//}
		eatcomment(inp);
	}
	return markermap;
}

LinkageMap read_eval_pos_file(const string& filename)
{
	LinkageMap positions;
	string line;
	ifstream inp;
	OpenFile(inp,filename);
	int line_nr = 0;
	while (getline(inp,line))
	{
		line_nr++;
		if (line.empty()) continue;
		istringstream line_stream(line);
		int chr;
		double pos;
		line_stream >> chr >> pos;
		positions.push_back(Locus(chr,pos,EVAL_POS));
	}
	sort(positions.begin(),positions.end());
	return positions;
}

vector<double> make_rec_map_inbred_lines(const LinkageMap& markermap)
{
	vector<double> r = make_rec_map(markermap);
	int M = r.size();
	for (int i=0;i<M;i++)
		r[i] = 2*r[i]/(1+2*r[i]);
	return r;
}

/*
// here we assume that the loci in the markermap are ordered!
int NrChromosomes(const LinkageMap& markermap)
{
	return markermap.back().GetChr() + 1;
}

vector<bool> get_selected_chr(const Args& Argu, const LinkageMap& markermap)
{
	int chr;
	const int Nchr = NrChromosomes(markermap);
	vector<bool> sel_chr(Nchr,true);
	if (Argu.GetOption("chr",chr))
	{
		if (chr < 1 || chr > Nchr)
			throw mblib_error("Error in -chr option (wrong chromosome number)");
		for (int i=0;i<Nchr;i++)
			sel_chr[i] = (i==chr-1) ? true : false;
	}
	return sel_chr;
}
*/

/*
LinkageMap add_small_dist_to_markermap(const LinkageMap& markermap,
									   const vector<bool>& sel_chr)
{
	const int M = markermap.size();
	const double eps = 1.0e-6;
	const double delta = 0.001;
	LinkageMap new_markermap;
	for (int m=0;m<M;m++)
	{
		Locus loc = markermap[m];
		if (sel_chr[loc.GetChr()])
		{
			if (!new_markermap.empty())
			{
				if (loc.GetChr() == new_markermap.back().GetChr())
				{
					double dist = loc.GetPosition() - new_markermap.back().GetPosition();
					if (dist < eps)
					{
						string n = loc.GetName();
						double p = new_markermap.back().GetPosition();
						int c = loc.GetChr();
						loc = Locus(c,p+delta,n);
					}
				}
			}
			new_markermap.push_back(loc);
		}
	}
	return new_markermap;
}
*/

