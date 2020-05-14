// Martin Boer, Biometris
#include <list>
#include <set>
#include <map>
#include <numeric>
#include <cmath>
#include <fstream>
#include <sstream>
#include "read_loc.h"
#include "read_map.h"

using namespace mbl;
using namespace std;

/*
int read_loc_file(vector<string>& ID,
				  vector<string>& marker_names,
				  matrix<score>& geno_obs,
				  const string& locfile)
{
	ifstream inp;
	OpenFile(inp,locfile);

	string line;
	getline(inp,line);

	// 12 jan 2007, here I removed the skipf() statement.
	std::istringstream line_stream(line);
	string tmp;
	line_stream >> tmp;
	int M;
	for (int cnt=0;;cnt++)
	{
		string name;
		line_stream >> name;
		if (!line_stream)
		{
			M = cnt;
			break;
		}
		tolower(name);
		marker_names.push_back(name);
	}
	while (getline(inp,line))
	{
		string id_ind;
		if (line.empty()) continue;
		istringstream line_stream(line);
		line_stream >> id_ind;
		ID.push_back(id_ind);
		vector<score> obs_ind;
		for (int m=0;m<M;m++)
		{
			score sc = read_score(line_stream);
			obs_ind.push_back(sc);
		}
		geno_obs.push_back(obs_ind);
	}
	return 0;
}
*/

// convert old format to new format,
// e.g. a score: 1,2 --> (1,2)
int convert_loc_file()
{
	ifstream inp;
	ofstream outp;
	OpenFile(inp,"Martin_geno.dat");
	OpenFile(outp,"Martin_geno_new.dat");
	outp.setf(ios_base::left,ios_base::adjustfield);

	string line;
	getline(inp,line);

	// 12 jan 2007, here I removed the skipf() statement.
	std::istringstream line_stream(line);
	string tmp;
	line_stream >> tmp;
	outp << setw(15) << tmp;
	int M;
	for (int cnt=0;;cnt++)
	{
		string name;
		line_stream >> name;
		if (!line_stream)
		{
			M = cnt;
			break;
		}
		outp << setw(15) << name;
	}
	outp << endl;
	while (getline(inp,line))
	{
		string id_ind;
		if (line.empty()) continue;
		istringstream line_stream(line);
		line_stream >> id_ind;
		outp << setw(15) << id_ind;
		for (int m=0;m<M;m++)
		{
			string score,pr;
			line_stream >> score;
			if (score == "*")
				pr = "  " + score + "  ";
			else
				pr = "(" + score + ")";
			outp << setw(15) << pr;
		}
		outp << endl;
	}
	return 0;
}



