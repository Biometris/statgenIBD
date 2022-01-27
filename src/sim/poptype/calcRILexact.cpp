#include "calcRILexact.h"
#include <bitset>
#include <map>

using namespace ibd;
using namespace std;

int ibd::countbits(unsigned int x)
{
	int cnt=0;
	for (;x !=0; x>>=1)
		if (x & 01)
			cnt++;
	return cnt;
}

bool ibd::shift_right(unsigned int& x)
{
	bool bit = (x & 01);
	x >>= 1;
	return bit;
}


Genotype get_genotype(InhVec u, int Nind)
{
	vector<char> c_new(2);
	vector<char> c = init_vector<char>('A','B');
	
	for (int i=0;i<Nind;i++)
	{
		c_new[0] = c[shift_right(u)];
		c_new[1] = c[shift_right(u)];
		c = c_new;
	}
	return Genotype(c[0],c[1]);
}

int nrecomb(const InhVec& v1, const InhVec& v2)
{
	return countbits(v1^v2);
}


Group find_group(InhVec u, int Nind)
{
	Group x = 0;

	vector<char> c_new(2);
	vector<char> c = init_vector<char>('A','B');
	
	for (int i=0;i<Nind;i++)
	{
		int k1 = shift_right(u);
		int k2 = shift_right(u);
		if (k1 != k2)
			x += pow2(i+1);

		c_new[0] = c[k1];
		c_new[1] = c[k2];
		c = c_new;
	}
	if (c[0] != 'A' && c[1] != 'A')
		x += 1;
	return x;
}

InhVec find_inhvec(Group group, int Nind)
{
	unsigned int M = pow2(Nind*2);
	for (InhVec u=0;u<M;u++)
		if (find_group(u,Nind) == group)
			return u;
	return 0;
}


class TransMatrixRI
{
public:
	TransMatrixRI(int x) : Ngen(x-1),Ngrp(pow2(x)-1)
	{
		A = matrix3D<int>(Ngrp,Ngrp,2*Ngen+1,0);
		unsigned int M = pow2(2*Ngen);
		
		for (Group i=0;i< Ngrp;i++)
		{
			InhVec u1 = find_inhvec(i,Ngen);
			for (InhVec u2=0;u2<M;u2++)
			{
				Group j = find_group(u2,Ngen);
				int nrec = nrecomb(u1,u2);
				A[i][j][nrec]++;
			}			
		}
	}

	matrix<double> operator()(double r) const
	{
		int K = 2*Ngen;
		double s = 1-r;
		vector<double> table(K+1); 
		for (int k=0;k<=K;k++)
			table[k] = pow(r,k)*pow(s,K-k);

		matrix<double> T(Ngrp,Ngrp);
		for (int i=0;i<Ngrp;i++)
		{
			for (int j=0;j<Ngrp;j++)
			{
				double sum = 0.0;
				for (int k=0;k<=K;k++)
					sum += A[i][j][k]*table[k];
				T[i][j] = sum;
			}
		}
		return T;
	}

private:
	int Ngrp, Ngen;
	matrix3D<int> A;
};


bool test_score(const Genotype& g, const ObsGeno& obs)
{
	for (ObsGeno::const_iterator it=obs.begin();it!=obs.end();it++)
		if (g == *it)
			return true;
	return false;
}

ObsInhVec obs_inh_vectors_aux(const ObsGeno& g, int Nind)
{
	ObsInhVec result;
	unsigned int M = pow2(2*Nind);
	for (InhVec k=0;k<M;k++)
	{
		if (test_score(get_genotype(k,Nind),g))
			result.push_back(k);
	}
	return result;
}

vector<ObsInhVec> ibd::obs_inh_vectors(const vector<ObsGeno>& g, int Nind)
{
	int nloc = g.size();
	vector<ObsInhVec> result(nloc);
	for (int m=0;m<nloc;m++)
		result[m] = obs_inh_vectors_aux(g[m],Nind);
	return result;
}


double transprob(const InhVec& v1, const InhVec& v2, int M, double r)
{
	const int k = nrecomb(v1,v2);
	double s = 1.0 - r;
	return pow(r,k)*pow(s,M-k);
}

matrix<double> transmatrix(int M, double r)
{
	unsigned int dim = pow2(M);
	matrix<double> T(dim,dim);
	for (InhVec i=0;i<dim;i++)
		for (InhVec j=0;j<dim;j++)
			T[i][j] = transprob(i,j,M,r);
	return T;
}

matrix<double> ibd::calc_prob_left(const vector<ObsInhVec>& obs_inh_vec, int M,
							  const vector<double>& r)
{
	const int nloc = obs_inh_vec.size();
	const int dim = pow2(M);
	matrix<double> L(nloc,dim,0.0);

	ObsInhVec ndx = obs_inh_vec[0];
	for (unsigned int k=0;k<ndx.size();k++)
		L[0][ndx[k]] = transprob(0,ndx[k],M,0.5);
	make_conditional(L[0]);
	for (int loc=1;loc<nloc;loc++)
	{
		ndx = obs_inh_vec[loc];
		for (unsigned int k=0;k<ndx.size();k++)
			for (InhVec q=0;q<(unsigned)dim;q++)
				L[loc][ndx[k]] += L[loc-1][q]*transprob(q,ndx[k],M,r[loc-1]);
		make_conditional(L[loc]);
	}
	return L;
}

matrix<double> ibd::calc_prob_right(const vector<ObsInhVec>& obs_inh_vec, int M,
							  const vector<double>& r)
{
	const int nloc = obs_inh_vec.size();
	const int dim = pow2(M);
	matrix<double> R(nloc,dim,0.0);

	ObsInhVec ndx = obs_inh_vec[nloc-1];
	for (unsigned int k=0;k<ndx.size();k++)
		R[nloc-1][ndx[k]] = transprob(ndx[k],0,M,0.5);
	make_conditional(R[nloc-1]);
	for (int loc=nloc-2;loc>=0;loc--)
	{
		ndx = obs_inh_vec[loc];
		for (unsigned int k=0;k<ndx.size();k++)
			for (InhVec q=0;q<(unsigned)dim;q++)
				R[loc][ndx[k]] += R[loc+1][q]*transprob(ndx[k],q,M,r[loc]);
		make_conditional(R[loc]);
	}
	return R;
}


vector<double> calc_QTL_prob_ind(const ObsGeno& U,
							 const vector<double>& L,
							 const vector<double>& R,
							 const matrix<double>& T1,
							 const matrix<double>& T2,
							 int M)
{
	const int Nind = M/2;
	map<Genotype,int> ndx;
	for (unsigned int q=0;q<U.size();q++)
		ndx[U[q]] = q;

	unsigned int dim = pow2(M);
	vector<double> left_prob(dim,0.0);
	vector<double> right_prob(dim,0.0);
	for (InhVec q1=0;q1<dim;q1++)
	{
		for (InhVec q2=0;q2<dim;q2++)
		{
			left_prob[q1] += L[q2]*T1[q2][q1];
			right_prob[q1] += T2[q1][q2]*R[q2];
		}
	}

	vector<double> prob(U.size(),0.0);
	for (InhVec q=0;q<dim;q++)
	{
		Genotype g = get_genotype(q,Nind);
		prob[ndx[g]] += left_prob[q]*right_prob[q];
	}
	make_conditional(prob);
	return prob;
}

matrix<double> ibd::calc_QTL_prob(const ObsGeno& U,
							 const LinkageMap& markermap,
							 const matrix3D<double>& Tr_l,
							 const matrix3D<double>& Tr_r,
							 int nind, int M,
							 const Locus& QTLpos)
{
	const int dimU = U.size();
	matrix<double> P(nind,dimU);
	const int left = pos_qtl(markermap,QTLpos);
	const int right = left + 1;
	
	double r_left = recomb(markermap[left],QTLpos);
	double r_right = recomb(QTLpos,markermap[right]);

	matrix<double> T1 = transmatrix(M,r_left);
	matrix<double> T2 = transmatrix(M,r_right);

	for (int i=0;i<nind;i++)
	{
		const vector<double>& L = Tr_l[i][left];
		const vector<double>& R = Tr_r[i][right];
		P[i] = calc_QTL_prob_ind(U,L,R,T1,T2,M);
	}

	return P;
}


// need to check!! Problem with conversion from Visual 2005 to 2010.
// It seems that this function is only used for test output.
vector<bool> conv_int_vecbl(const unsigned int& x, int n)
{
	unsigned long long z = x;
	const bitset<32> y = z;
	vector<bool> result(n);
	for (int i=0;i<n;i++)
		result[i] = y[n-i-1];
	return result;
}

vector<Genotype> make_genotype(InhVec u, int Nind)
{
	vector<Genotype> G(Nind);

	vector<char> c_new(2);
	vector<char> c = init_vector<char>('A','B');
	
	for (int i=0;i<Nind;i++)
	{
		c_new[0] = c[shift_right(u)];
		c_new[1] = c[shift_right(u)];
		c = c_new;
		G[Nind-i-1] = Genotype(c[0],c[1]);
	}
	return G;
}

void print_recomb(int Nind,double r)
{
	int Ngam = 2*Nind;
	int M = pow2(Ngam);
	for (InhVec i=0;i<M;i++)
	{
		for (InhVec j=0;j<M;j++)
		{
			cout << conv_int_vecbl(i,Ngam) << "  " 
				 << conv_int_vecbl(j,Ngam) << "  " 
			     << setw(3) << make_genotype(i,Nind) 
				 << setw(3) << make_genotype(j,Nind)
				 << setw(6) << nrecomb(i,j) 
				 << setw(12) << transprob(i,j,Ngam,r) << endl;
		}
	}
}

int ibd::test_calc_prob_RIL()
{
	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(3) << endl;
	TransMatrixRI T(3);
	cout << setw(6) << T(0.10) << endl;

	return 0;

	for (InhVec k=0;k<pow2(2*3);k++)
	{
		cout << conv_int_vecbl(k,2*3) << "     "  
			 << setw(3) << find_group(k,3) << endl;
	}

	return 0;


	cout.setf(ios::fixed, ios::floatfield);

	double r1 = 0.25;
	double r2 = 0.10;
	double r12 = r1*(1-r2) + r2*(1-r1);

	matrix<double> H1 = transmatrix(6,r1);
	matrix<double> H2 = transmatrix(6,r2);
	matrix<double> H12 = transmatrix(6,r12);
	matrix<double> Chk = H1*H2 - H12;
	cout << "DIM: " << Chk.NrRows() << endl;
	for (unsigned int i=0;i!=Chk.NrRows();i++)
		for (unsigned int j=0;j!=Chk.NrRows();j++)
			if (Chk[i][j] > 1.0e-8)
				cout << setw(3) << i << setw(3) << j << endl;

	
	const Genotype A('A','A');
	const Genotype H('A','B');
	const Genotype B('B','B');
	const ObsGeno C(B,H);
	const ObsGeno D(A,H);
	const ObsGeno U(A,H,B);

	map<Genotype,int> cnt;
	const int Ngen = 2;
	const int Ngam = 2*Ngen;
	unsigned int M = pow2(Ngam);
	for (InhVec k=0;k<M;k++)
	{
		vector<Genotype> geno = make_genotype(k,Ngen);
		cout << conv_int_vecbl(k,Ngam) << "     "  << setw(3) << geno;
		Genotype g = geno[0];
		cnt[g]++;
		if (test_score(g,C))
			cout << " * ";
		cout << endl;
	}
	cout << "# poss.: " << M << endl;
	for (map<Genotype,int>::const_iterator it=cnt.begin();it!=cnt.end();it++)
		cout << setw(3) << it->first << setw(12) << it->second << endl;

	ObsInhVec obs = obs_inh_vectors_aux(C,Ngen);
	vector<ObsGeno> g(2,C);
	vector<ObsInhVec> tst = obs_inh_vectors(g,Ngen);

	cout << "Obs.size() = " << obs.size() << endl;
	for (unsigned int i=0;i<obs.size();i++)
		cout << setw(3) << i << "  " 
		     << conv_int_vecbl(obs[i],Ngam) << " "
			 << setw(3) << make_genotype(obs[i],Ngen) << " "
			 << setw(3) << make_genotype(tst[0][i],Ngen) << " "
			 << setw(3) << make_genotype(tst[0][i],Ngen) << endl;

	print_recomb(Ngen,0.10);
	return 0;
}

/*
InhVec conc(const InhVec& x,bool p1,bool p2)
{
	InhVec y = x;
    y.push_back(p1);
	y.push_back(p2);
	return y;
}

void bin_tree(int Ngen, const InhVec& x,vector<InhVec>& C)
{
	if (Ngen == 0) 
		C.push_back(x);
	else
	{
		bin_tree(Ngen-1,conc(x,0,0),C);
		bin_tree(Ngen-1,conc(x,0,1),C);
		bin_tree(Ngen-1,conc(x,1,0),C);
		bin_tree(Ngen-1,conc(x,1,1),C);
	}
}

unsigned int conv_vecbl_int(const vector<bool>& x)
{
	unsigned int result = 0;
	unsigned int a=1;
	for (int i=x.size()-1;i>=0;i--)
	{
		if (x[i])
			result += a;
		a*=2;
	}
	return result;
}

bool getbit(unsigned int x, int p)
{
	return (x >> p) & 01;
}

*/

