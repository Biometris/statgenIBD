// Martin Boer, Biometris
#include <iostream>
#include <sstream>
#include <utility> 
#include <fstream>
#include <algorithm>
#include <locale>
#include "misc.h"
#include "mblexcept.h"
#include "convert.h"
#include "markerscore.h"

using namespace std;
using namespace mbl;

const bool scores_in_parentheses = true; 

score::score(int a, int b) : pair<int,int>(a,b) 
{ 
	if (a < b) 
		std::swap(this->first,this->second);
}

bool score::homozygous() const
{
	return !(this->first == U_haplo_sc || this->second == U_haplo_sc || 
			 this->first != this->second);
}

string score::print_string() const
{
	string f="*";
	string s="*";
	if (this->first != U_haplo_sc)
		f = stringify(this->first);
	if (this->second != U_haplo_sc)
		s = stringify(this->second);
	return "(" + f + "," + s + ")";
}

std::ostream& operator<<(std::ostream& outp, const score& score)
{
	string str = score.print_string();
	int nr_whitespaces = outp.width() - str.length();
	outp.width(0);
	if (nr_whitespaces > 0)
		outp << string(nr_whitespaces,' ');
	outp << str;
	return outp;
}


// reads a character from input and checks whether this character is equal to c
void check_char(istream& s, char c)
{
	char x;
	s >> x;
	if (x != c)
		throw mblib_error("missing " + string(1,c));
}

int read_allele(istream& s)
{
	char c;
	s >> c;
	if (c == '-')
		return U_haplo_sc;	
	s.putback(c);

	char x;
	if (s >> x)
	{
		//if (x < 0)
		//	throw mblib_error("allele value negative!");
		//else
		return int(x);
	}
	
	throw mblib_error("error while reading allele");
	return 0; // dummy
}

// temp. added, april 2013
// assume scores 'a', 'b' and '-':
//score read_score(istream& s)
//{
//	char c;
//	s >> c;
//	if (c == '-')
//		return Uscore;
//	if (c == 'a')
//		return score(1,1);
//	if (c == 'b')
//		return score(2,2);
//}

score read_score(istream& line_stream, char delimit)
{
  string field;
  getline(line_stream,field,delimit);
  field.erase(remove_if(field.begin(), field.end(),::isspace), field.end()); 
    
  istringstream s(field);
  if (field == "-") return Uscore;

  string::size_type pos = field.find('/');
	if (pos == field.npos)     // score is of the form "a"
	{
     int a = read_allele(s);
     return score(a,a);
	}

  int a = read_allele(s);
	check_char(s,'/');
	int b = read_allele(s);

	return score(a,b);
}

bool check_score(const OrdGeno& g, const score& sc)
{
	if (sc == Uscore) 
		return true;
	if (sc.second == U_haplo_sc)
		return (sc.first == g.first || sc.first == g.second);
	return (((g.first == sc.first)&&(g.second == sc.second))||
			((g.first == sc.second)&&(g.second == sc.first)));
}

vector<string> read_tab_delimited_line(istream& inp)
{
  vector<string> result;
	string line,str;
	getline(inp,line);
  istringstream line_stream(line);
	while (getline(line_stream,str,'\t'))
	{
		result.push_back(str);
	}
	return result;
}

int read_flapjackfile(vector<string>& geno, vector<string>& markers,matrix<score>& scores,
      const string filename)
{
  ifstream inp;
  OpenFile(inp,filename);  

  string str;
  getline(inp,str,'\t');
	markers = read_tab_delimited_line(inp); 
	
  const int M = markers.size();
	string line;
	while (getline(inp,line))
	{
		if (line.empty()) break;
    istringstream line_stream(line);
    string g;
    getline(line_stream,g,'\t');
    geno.push_back(g);
    
    vector<score> ind_scores;
    for (int m=0;m<M;m++)
    {
       score sc = read_score(line_stream,'\t');
       ind_scores.push_back(sc);
    }
    scores.push_back(ind_scores);
	}
  return 0;
}

/*
int test_read()
{
	//				    A      H     B     C     D   U
	istringstream inp("(1,1) (1,2) (2,2) (2,*) (1,*) *");
	const int N = 6;
	cout << "ordered geno:  [1,1] [1,2] [2,1] [2,2] " << endl;

	for (int i=0;i<N;i++)
	{
		score sc = read_score(inp);
		cout << setw(4) << i << setw(8) << sc;
		for (int a1=1;a1<=2;a1++)
		{
			for (int a2=1;a2<=2;a2++)
			{
				OrdGeno g(a1,a2);
				cout << setw(6) << check_score(g,sc);
			}
		}
		cout << endl;
	}
	return 0;
}
*/

vector<int> haplo_score(const vector<score>& obsgeno)
{
	int M = obsgeno.size();
	vector<int> haplo(M);
	for (int j=0;j<M;j++)
	{
		score g = obsgeno[j];
		haplo[j] = g.homozygous() ? g.first : U_haplo_sc;
	}
	return haplo;
}

matrix<int> 
gen_haplo_score_ped(vector<IndProp>& sel_pop,
				 const matrix<score>& obsgeno, 
				 const vector<IndProp>& pop)
{
	const int N = obsgeno.NrRows();
	matrix<int> result;
	for (int i=0;i<N;i++)
	{
		const IndProp& ind = pop[i];
		if (!ind.IsMemberFamily())
		{
			vector<int> haplo = haplo_score(obsgeno[i]);
			result.push_back(haplo);
			sel_pop.push_back(ind);
		}
	}
	return result;
}

// haplo scores for families + parents families
matrix<score>
select_families(vector<IndProp>& sel_pop,
				const matrix<score>& obsgeno, 
				const vector<IndProp>& pop, 
				const set<string>& Pfam)
{
	const int N = obsgeno.NrRows();
	matrix<score> result;
	for (int i=0;i<N;i++)
	{
		IndProp ind = pop[i];
		bool parent_fam = (Pfam.find(ind.GetID()) != Pfam.end());
		if (ind.IsMemberFamily() || parent_fam)
		{
			string ID   = ind.GetID();
			string type = ind.GetType();
			string fam  = ind.GetFam();
			string P1   = ind.GetP1();
			string P2   = ind.GetP2();
			if (parent_fam)
			{
				type = "INBFND";
			}
			result.push_back(obsgeno[i]);
			sel_pop.push_back(IndProp(ID,fam,type,P1,P2));
		}
	}
	return result;
}

vector<bool> detect_high_density_ind(const matrix<int>& geno, int min)
{
	int N = geno.NrRows();
	int M = geno.NrCols();
	vector<bool> hd(N);
	for (int i=0;i<N;i++)
	{
		int cnt = 0;
		for (int m=0;m<M;m++)
			if (geno[i][m] != -1)
				cnt++;
		hd[i] = (cnt >= min);		
	}
	return hd;
}

// npar: number of parents
map<score,int> ndx_score(int npar)
{
	int k=0;
	map<score,int> ndx;
	for (int i=0;i<npar;i++)
		ndx[score(i,i)] = k++;
	if (npar == 2)
	{
		ndx[score(0,1)] = k++;
	}
	else if (npar == 3)
	{
		for (int i=0;i<2;i++)
			ndx[score(i,2)] = k++;
	}
	else if (npar == 4)
	{
		for (int i=0;i<2;i++)
			for (int j=2;j<4;j++)
				ndx[score(i,j)] = k++;
	}
	else
		throw mblib_error("npar > 4 !");
	return ndx;
}

int tst_ndx_score()
{
	ofstream outp;
	OpenFile(outp,"tmp.dat");
	for (int npar = 2; npar < 5; npar++)
	{
		outp << "npar: " << npar << endl;
		map<score,int> tst = ndx_score(npar);
		for (map<score,int>::const_iterator it=tst.begin();it!=tst.end();it++)
		{	
			outp << it->first << '\t' << it->second << endl;
		}
		outp << endl;
	}
	return 0;
}



/*
// haplo scores for families + parents families
matrix<int>
gen_haplo_score_fam(vector<IndProp>& sel_pop,
					const matrix<score>& obsgeno, 
				    const vector<IndProp>& pop, 
				    const set<string>& Pfam)
{
	const int N = obsgeno.NrRows();
	matrix<int> result;
	for (int i=0;i<N;i++)
	{
		IndProp ind = pop[i];
		bool parent_fam = (Pfam.find(ind.GetID()) != Pfam.end());
		if (ind.IsMemberFamily() || parent_fam)
		{
			string ID   = ind.GetID();
			string type = ind.GetType();
			string P1   = ind.GetP1();
			string P2   = ind.GetP2();
			vector<int> haplo = haplo_score(obsgeno[i]);
			if (parent_fam)
			{
				type = "INBFND";
				// P1 = P2 = "0"; // founder
			}
			result.push_back(haplo);
			sel_pop.push_back(IndProp(ID,type,P1,P2));
		}
	}
	return result;
}
*/

/*
score read_score(istream& s)
{
	score result;
	char c;
	s >> c;
	if (c == '*') return Uscore;

	s.putback(c);
	s >> result.first >> c >> result.second;
	if (result.first < 0 || result.second < 0 || c != ',')
	   throw mblib_error("Error while reading score!!: ");
	return result;
}
*/

