#include <algorithm>
#include <iostream>
#include <fstream>

#include "convert.h"
#include "ibdexcept.h"
#include "util_genetics.h"
#include "misc.h"

using namespace ibd;
using namespace std;

IndProp::IndProp(std::string id,
                 std::string fam,
                 std::string type,
                 std::string par1,
                 std::string par2)
	: ID(id), FAM(fam), TYPE(type), P1(par1), P2(par2) { }

// simple function object (functor), used by findID function.
class eqID
{
public:
	eqID(const string& id) : ID(id) {}
	bool operator()(const IndProp& ind) { return ind.GetID() == ID; }
private:
	string ID;
};

bool findID(const vector<IndProp>& pop, const string& ID)
{
	return (find_if(pop.begin(),pop.end(),eqID(ID)) != pop.end());
}

int ndxID(const vector<IndProp>& pop,
          const string& ID)
{
	return (find_if(pop.begin(),pop.end(),eqID(ID)) - pop.begin());
}

bool match(int& x,
           const string& str,
           const char * pat)
{
	string pattern(pat);
	int cnt = count(pattern.begin(),pattern.end(),'x');
	if (cnt != 1)
		throw ibd_error("error in match!");
	if (str.size() != pattern.size())
		return false;
	else
	{
		for (int i=0;i<(int)str.size();i++)
		{
			if (pattern[i] != 'x')
			{
				if (str[i]!=pattern[i])
					return false;
			}
			else
			{
				string tmp = string(1,str[i]);
				x = convertTo<int>(tmp);
			}
		}
	}
	return true;
}

// function used e.g. by BCxSyDH population
bool match(int& x,
           int& y,
           const string& str,
           const char * pat)
{
	string pattern(pat);
	int cnt_x = count(pattern.begin(),pattern.end(),'x');
	int cnt_y = count(pattern.begin(),pattern.end(),'y');
	if (cnt_x != 1)
		throw ibd_error("error in match!");
	if (cnt_y != 1)
		throw ibd_error("error in match!");

	if (str.size() != pattern.size())
		return false;
	else
	{
		for (int i=0;i<(int)str.size();i++)
		{
			if (pattern[i] == 'x')
			{
				string tmp = string(1,str[i]);
				x = convertTo<int>(tmp);
			}
			else if (pattern[i] == 'y')
			{
				string tmp = string(1,str[i]);
				y = convertTo<int>(tmp);
			}
			else // pattern[i] != x && y
			{
				if (str[i]!=pattern[i])
					return false;
			}
		}
	}
	return true;
}

// use another name for this function!
bool correct_type(string type)
{
	toupper(type);
	return (type == "INBFND" || type == "INBPAR" || type == "RIL" || type == "HYBRID");
}

vector<IndProp> read_ped_file(const string& filename)
{
	vector<IndProp> pop;
	string line;
	std::ifstream inp;
	OpenFile(inp,filename);
	int line_nr = 0;
	while (getline(inp,line))
	{
		line_nr++;
		if (line.empty()) continue;
		istringstream line_stream(line);
		string ID,type,fam,P1,P2;
		if (fam_column_pedigree_file)
			line_stream >> ID >> fam >> type >> P1 >> P2;
		else
		{
			line_stream >> ID >> type >> P1 >> P2;
			if (correct_type(type))
				fam = "*";
			else
				fam = "ID_FAM";
		}
		if (findID(pop,ID))
			throw ibd_error("ID " + ID + " not unique");
		if (fam == "*" && !correct_type(type))
			throw ibd_error("type " + type + " not defined");

		if (type != "INBFND" && type != "INBPAR")
		{
			if (!findID(pop,P1))
				throw ibd_error("Parent " + P1 + " not defined");
			if (!findID(pop,P2))
				throw ibd_error("Parent " + P2 + " not defined");
		}
		pop.push_back(IndProp(ID,fam,type,P1,P2));
	}
	return pop;
}

