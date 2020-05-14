#include <algorithm>
#include <numeric>
#include <fstream>
#include "eigen.h"
#include "linalg.h"
#include "misc.h"

#include "util.h"
//#include "Random.h"

using namespace std;
using namespace mbl;

namespace
{
    inline void rot(matrix<double>& a, double s, double tau, int i, int j, int k, int l)
	{
		double g = a[i][j];
		double h = a[k][l];
		a[i][j] = g - s*(h+g*tau);
		a[k][l] = h + s*(g-h*tau);
	}
}

void mbl::Jacobi(matrix<double> a, vector<double>& d, matrix<double>& v)
{
	int j,ip,iq;
	double tresh,theta,tau,t,sm,s,h,g,c;

	const int n = a.NrCols();
	vector<double> b(n),z(n);
	v = identity_matrix(n);
	d = vector<double>(n);
	for (ip=0;ip<n;ip++)
	{
		b[ip] = d[ip] = a[ip][ip];
        z[ip] = 0.0;
	}
	for (int i=1;i<=50;i++)
	{
		sm = 0.0;
		for (ip=0;ip<n-1;ip++)
		{
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) return;
		if (i < 4)
			tresh = 0.2*sm/(n*n);
		else
			tresh = 0.0;
		for (ip=0;ip<n-1;ip++)
		{
			for (iq=ip+1;iq<n;iq++)
			{
				g = 100.0*fabs(a[ip][iq]);
				if (i>4 && (fabs(d[ip])+g) == fabs(d[ip]) && (fabs(d[iq])+g) == fabs(d[iq]))
					a[ip][iq] = 0.0;
				else if (fabs(a[ip][iq]) > tresh)
				{
					h = d[iq]-d[ip];
					if ((fabs(h)+g) == fabs(h))
						t = a[ip][iq]/h;
					else
					{
						theta = 0.5*h/a[ip][iq];
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c = 1.0/sqrt(1+t*t);
					s = t*c;
					tau = s/(1.0+c);
					h = t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
					for (j=0;j<ip;j++)
						rot(a,s,tau,j,ip,j,iq);
					for (j=ip+1;j<iq;j++)
						rot(a,s,tau,ip,j,j,iq);
					for (j=iq+1;j<n;j++)
						rot(a,s,tau,ip,j,iq,j);
					for (j=0;j<n;j++)
						rot(v,s,tau,j,ip,j,iq);
				}
			}
		}
		for (ip=0;ip<n;ip++)
		{
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	throw mblib_error("Too many iterations in routine Jacobi");
	return;
}

// calculation of eigenvalues for a SYMMETRIC matrix A
vector<double> mbl::eigenvals(matrix<double>& A)
{
	vector<double> d;
	matrix<double> v;
	Jacobi(A,d,v);
	sort(d.begin(),d.end());
	return d;
}

// square root of a SYMMETRIC MATRIX M
// sqrt_M = u*sqrt_d*v,  v = u^T, v*u = I
// sqrt_M * sqrt_M = (u*sqrt_d*v)*(u*sqrt_d*v) = u*d*v = M
matrix<double> mbl::square_root(const matrix<double>& M)
{
	int dim = M.NrCols();
	vector<double> d; // eigenvalues
	matrix<double> u; // normalized eigenvectors in the columns
    Jacobi(M, d, u);  // M = u*d*u^T = u*d*v

	vector<double> sqrt_d(dim);
	for (int k=0;k<dim;k++)
		sqrt_d[k] = sqrt(d[k]);

	matrix<double> sqrt_M(dim,dim);
	for (int i=0;i<dim;i++)
	{
		for (int j=0;j<dim;j++)
		{
			double sum = 0.0;
			for (int k=0;k<dim;k++)
				sum += sqrt_d[k]*u[i][k]*u[j][k];
			sqrt_M[i][j] = sum;
		}
	}
	return sqrt_M;
}



// deze functies eventueel nog vervangen door standaardfuncties
inline double sign(double a, double b) { return b > 0.0 ? fabs(a) : -fabs(a); }
inline int max(int a, int b) { return a > b ? a : b; }
inline void error(const char *msg)
{
	cerr << "Error: " << msg << endl;
	exit(1);
}

void balanc(matrix<double>& a)
{
	const double RADIX = numeric_limits<double>::radix;
	int i,j,last=0;
	double s,r,g,f,c,sqrdx;

	int n=a.NrRows();
	sqrdx = RADIX*RADIX;
	while (last == 0)
	{
		last = 1;
		for (i=0;i<n;i++)
		{
			r = c = 0.0;
			for (j=0;j<n;j++)
			{
				if (j!=i)
				{
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			}
			if (c != 0.0 && r != 0.0)
			{
				g = r/RADIX;
				f = 1.0;
				s = c + r;
				while (c < g)
				{
					f *= RADIX;
					c *= sqrdx;
				}
				g = r*RADIX;
				while (c > g)
				{
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c+r)/f < 0.95*s)
				{
					last = 0;
					g = 1.0/f;
					for (j=0;j<n;j++)
						a[i][j] *= g;
					for (j=0;j<n;j++)
						a[j][i] *= f;

				}
			}
		}
	}
}


void elmhes(matrix<double>& a)
{
	int i,j,m;
	double y,x;

	int n=a.NrRows();
	for (m=1;m<n-1;m++)
	{
		x = 0.0;
		i = m;
		for (j=m;j<n;j++)
		{
			if (fabs(a[j][m-1]) > fabs(x))
			{
				x = a[j][m-1];
				i = j;
			}
		}
		if (i != m)
		{
			for (j=m-1;j<n;j++)
				swap(a[i][j],a[m][j]);
			for (j=0;j<n;j++)
				swap(a[j][i],a[j][m]);
		}
		if (x != 0.0)
		{
			for (i=m+1;i<n;i++)
			{
				if ((y=a[i][m-1]) != 0.0)
				{
					y /= x;
					a[i][m-1] = y;
					for (j=m;j<n;j++)
						a[i][j] -= y*a[m][j];
					for (j=0;j<n;j++)
						a[j][m] += y*a[j][i];
				}
			}
		}
	}
}

vector< complex<double> > hqr(matrix<double>& a)
{
	int nn,m,l,k,j,its,i,mmin;
	double z,y,x,w,v,u,t,s,r,q,p,anorm;

	int n = a.NrRows();
	vector< complex<double> > wri(n);

	anorm = 0.0;
	for (i=0;i<n;i++)
		for (j=max(i-1,0);j<n;j++)
			anorm += fabs(a[i][j]);
	nn = n-1;
	t = 0.0;
	while (nn >= 0)
	{
		its = 0;
		do
		{
			for (l=nn;l>0;l--)
			{
				s = fabs(a[l-1][l-1]) + fabs(a[l][l]);
				if (s == 0.0) s = anorm;
				if (fabs(a[l][l-1]) + s == s)
				{
					a[l][l-1] = 0.0;
					break;
				}
			}
			x = a[nn][nn];
			if (l == nn)
			{
				wri[nn--] = x+t;
			}
			else
			{
				y = a[nn-1][nn-1];
				w = a[nn][nn-1]*a[nn-1][nn];
				if (l == nn-1)
				{
					p = 0.5*(y-x);
					q = p*p +w;
					z = sqrt(fabs(q));
					x += t;
					if (q >= 0.0)
					{
						z = p + sign(z,p);
						wri[nn-1] = wri[nn] = x + z;
						if (z != 0.0) wri[nn] = x-w/z;
					}
					else
					{
						wri[nn] = complex<double>(x+p,z);
						wri[nn-1] = conj(wri[nn]);
					}
					nn -=2;
				}
				else
				{
					if (its == 30) error("Too many iterations in hqr");
					if (its == 10 || its == 20)
					{
						t += x;
						for (i=0;i<nn+1;i++) a[i][i] -= x;
						s = fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
						y = x = 0.75*s;
						w = -0.4375*s*s;
					}
					++its;
					for (m=nn-2;m>=l;m--)
					{
						z = a[m][m];
						r = x-z;
						s = y-z;
						p = (r*s-w)/a[m+1][m]+a[m][m+1];
						q = a[m+1][m+1] - z - r - s;
						r = a[m+2][m+1];
						s = fabs(p) + fabs(q) + fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l) break;
						u = fabs(a[m][m-1])*(fabs(q)+fabs(r));
						v = fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
						if (u+v == v) break;
					}
					for (i=m;i<nn-1;i++)
					{
						a[i+2][i] = 0.0;
						if (i != m) a[i+2][i-1] = 0.0;
					}
					for (k=m;k<nn;k++)
					{
						if (k != m)
						{
							p = a[k][k-1];
							q = a[k+1][k-1];
							r = 0.0;
							if (k+1 != nn) r = a[k+2][k-1];
							if ((x=fabs(p)+fabs(q)+fabs(r)) != 0.0)
							{
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s = sign(sqrt(p*p+q*q+r*r),p)) != 0.0)
						{
							if (k == m)
							{
								if (l != m)
									a[k][k-1] = -a[k][k-1];
							}
							else
								a[k][k-1] = -s*x;
							p += s;
							x = p/s;
							y = q/s;
							z = r/s;
							q /= p;
							r /= p;
							for (j=k;j<nn+1;j++)
							{
								p = a[k][j] + q*a[k+1][j];
								if (k+1 != nn)
								{
									p+= r*a[k+2][j];
									a[k+2][j] -= p*z;
								}
								a[k+1][j] -= p*y;
								a[k][j] -= p*x;
							}
							mmin = nn < k+3 ? nn : k+3;
							for (i=l;i<mmin+1;i++)
							{
								p=x*a[i][k]+y*a[i][k+1];
								if (k != nn)
								{
									p += z*a[i][k+2];
									a[i][k+2] -= p*r;
								}
								a[i][k+1] -= p*q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l+1 < nn);
	}
	return wri;
}

struct order_eigenval
{
	order_eigenval(){}
	bool operator()(const complex<double>& x, const complex<double>& y)
	{
		if (real(x) == real(y))
			return (imag(x) < imag(y));
		return (real(x) > real(y));
	}
};

vector< complex<double> > mbl::eigenvalues(matrix<double> A)
{
	int i;
	balanc(A);
	elmhes(A);
	vector< complex<double> > result = hqr(A);
	const int n = A.NrCols();
	vector< complex<double> > res(n);
	for (i=0;i<n;i++)
		res[i] = result[i];

	sort(res.begin(),res.end(),order_eigenval());

	for (i=0;i<n;i++)
		result[i] = res[i];
	return result;
}


// Test of the calculation of the eigenvalues for a symmetric matrix. This
// example calculates the eigenvalues of the matrix Z^T P_0 Z for an unbalanced
// one-way ANOVA: K is the number of levels, J[i] is the number of
// observations per level and N is the total number of observations
// see also Crainiceanu & Ruppert, J. R. Statist. Soc. B (2004)
// Here we assume that Sigma = I.
void mbl::test_calc_eigenvals_sym()
{
	int i,j,k;
	vector<int> J = init_vector<int>(100,100,50);
	const int K = J.size();
	const int N = accumulate(J.begin(),J.end(),0);
	vector<double> res_var(N,1.0);
	/*
	for (i=0;i<N;i++)
		if (i < 100)
			res_var[i] = 1.0;
		else if (i < 200)
			res_var[i] = 1.0;
		else if (i < 300)
			res_var[i] = 8.0;
	*/
	matrix<double> R_inv(N,N,0.0);
	matrix<double> X(N,1,1.0);
	matrix<double> Z(N,K,0.0);
	for (i=0,k=0;k<K;k++)
	{
		for (j=0;j<J[k];j++)
		{
			Z[i][k] = 1.0;
			R_inv[i][i] = 1.0/res_var[i];
			i++;
		}
	}
	matrix<double> Xt = transpose(X);
	matrix<double> S0 = R_inv - R_inv*X*Inverse(Xt*R_inv*X)*Xt*R_inv;

	matrix<double> Q = transpose(Z)*S0*Z;

	vector<double> eigenvalues = eigenvals(Q);
	const int dim = eigenvalues.size();
	vector<double> det(dim);
	for (int d=0;d<dim;d++)
		det[d] = Determinant(Q-eigenvalues[d]*identity_matrix(dim));

	/*
	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(2);
	cout << "Eigenvals(Q):            " << setw(8) << eigenvalues << endl;
	cout << "Check -> det(Q-lambda*I):" << setw(8) << det << endl << endl;

	// calculation of effective dimension in two different ways
	double lambda = 100.0;
	matrix<double> X_all(N,K+1);
	for (i=0;i<N;i++)
	{
		X_all[i][0] = X[i][0];
		for (j=0;j<K;j++)
			X_all[i][j+1] = Z[i][j];
	}
	matrix<double> P = identity_matrix(K+1);
	P[0][0] = 0.0; // no penalty on main effect
	matrix<double> XtRinvX = transpose(X_all)*R_inv*X_all;
	matrix<double> H = Inverse(XtRinvX+lambda*P)*XtRinvX;
	cout << setprecision(6);

	double eff_dim1 = 0.0;
	for (i=0;i<H.size();i++)
	{
		cout << "H[i][i] " << H[i][i] << endl;
		eff_dim1 += H[i][i];
	}

	double eff_dim2 = 1.0;
	for (i=0;i<dim;i++)
	{
		cout << "mu_i/(lambda+mu_i)  " << eigenvalues[i]/(lambda+eigenvalues[i]) << endl;
		eff_dim2 += eigenvalues[i]/(lambda + eigenvalues[i]);
	}
	cout << "Eff Dim 1: "  << eff_dim1 << endl;
	cout << "Eff Dim 2: "  << eff_dim2 << endl;

	// test of calculation of asymptotic distribution of RLRT
	long int nsim = 100000;
	int ngridpoints = 200;
	vector<double> RLRT(nsim);

	matrix<double> A(ngridpoints,K);
	vector<double> b(ngridpoints,0.0);
	for (int s=0;s<ngridpoints;s++)
	{
		double d = exp(-12.0 + 24.0*s/(ngridpoints-1));
		for (int k=0;k<K;k++)
		{
			double dmu = d*eigenvalues[k];
			A[s][k] = dmu/(1+dmu);
			b[s] += log(1+dmu);
		}
	}
	int sim;
	for (sim=0;sim<nsim;sim++)
	{
		vector<double> w_sqr(K);
		for (int k=0;k<K;k++)
		    w_sqr[k] = sqr(randnormal());
		vector<double> RLRT_grid = A*w_sqr - b;
		double max = *max_element(RLRT_grid.begin(),RLRT_grid.end());
		RLRT[sim] = (max > 0.0) ? max : 0.0;
	}
	sort(RLRT.begin(),RLRT.end());
	ofstream outp;
	OpenFile(outp,"RLRT.dat");
	for (sim=0;sim<nsim;sim++)
		outp << setw(12) << (sim+1)/(1.0*nsim)
		     << setw(12) << RLRT[sim]
			 << setw(12) << sqrt(RLRT[sim]) << endl;
	cout << "Output written to file RLRT.dat" << endl << endl;
	*/
}

void mbl::test_calc_eigenvals_nonsym()
{
	cout << "Test of calculation of eigenvalues for non-sym matrix:" << endl;
	const int N=4;
	matrix<double> A(N,N);

	A[0][0] =  0.0; A[0][1] = -1.0; A[0][2] = 3.0; A[0][3] =  0.0;
	A[1][0] = -1.0; A[1][1] =  0.0; A[1][2] = 1.0; A[1][3] =  0.0;
	A[2][0] =  0.0; A[2][1] = -4.0; A[2][2] = 5.0; A[2][3] =  0.0;
	A[3][0] =  0.0; A[3][1] =  0.0; A[3][2] = 0.0; A[3][3] = -1.0;

	cout.setf(ios::fixed, ios::floatfield);		// set output format
	vector< complex<double> > eig = eigenvalues(A);
	cout << setw(3) << "nr" << setw(12) << "real" << setw(12) << "imag" << endl;
	for (int m=0;m<N;m++)
		cout << setw(3) << m+1 << setw(12) << real(eig[m]) << setw(12) << imag(eig[m]) << endl;
	cout << endl;
}

// added 26 januari 2009: Calculation of eigenvalues and eigenvectors.

// computes c=sqrt(a^2+b^2), without destructive underflow or overflow:
double pythag(double a, double b)
{
	double absa = fabs(a);
	double absb = fabs(b);
	if (absa > absb)
		return absa*sqrt(1.0+sqr(absb/absa));
	else
		return (absb = 0.0 ? 0.0 : absb*sqrt(1.0+sqr(absa/absb)));
}

void HouseHolderReduction(double ** a, double* d, double * e, int n)
{
	int i,j,l,k;
	double f,g,h,hh,scale;
	for (i=n-1;i>0;i--)
	{
		l=i-1;
		h = scale = 0.0;
		if (l > 0)
		{
			for (k=0;k<l+1;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i] = a[i][l];
			else
			{
				for (k=0;k<l+1;k++)
				{
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				f = a[i][l];
				g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale*g;
				h -= f*g;
				a[i][l] = f-g;
				f = 0.0;
				for (j=0;j<l+1;j++)
				{
					a[j][i]=a[i][j]/h;
					g = 0.0;
					for (k=0;k<j+1;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<l+1;k++)
						g += a[k][j]*a[i][k];
					e[j] = g/h;
					f += e[j]*a[i][j];
				}
				hh = f/(h+h);
				for (j=0;j<l+1;j++)
				{
					f = a[i][j];
					e[j] = g = e[j] - hh*f;
					for (k=0;k<j+1;k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		}
		else
			e[i]=a[i][l];
		d[i] = h;
	}
	d[0] = 0.0;
	e[0] = 0.0;
	for (i=0;i<n;i++)
	{
		l=i;
		if (d[i] != 0.0)
		{
			for (j=0;j<l;j++)
			{
				g = 0.0;
				for (k=0;k<l;k++)
					g += a[i][k]*a[k][j];
				for (k=0;k<l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i];
		a[i][i] = 1.0;
		for (j=0;j<l;j++)
			a[j][i] = a[i][j] = 0.0;
	}
}

void QL_algo(double * d, double * e, double ** z, int n)
{
	const int MAX_ITER = 250;
	int i,m,l,k,iter;
	double s,r,p,g,f,dd,c,b;

	for (i=1;i<n;i++)
		e[i-1] = e[i];
	e[n-1] = 0.0;
	for (l=0;l<n;l++)
	{
		iter = 0;
		do
		{
			for (m=l;m<n-1;m++)
			{
				dd = fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m]) + dd == dd)
					 break;

				//const double eps = numeric_limits<double>::epsilon();
				//double volatile temp = fabs(e[m]) + dd;
				//if (std::abs(e[m]) <= eps*dd)
				//	break;
				//if (temp == dd) break;
			}
			if (m != l)
			{
				if (iter++ == MAX_ITER)
					error("Too many iterations in QL_algo");
				g = (d[l+1]-d[l])/(2.0*e[l]);
				r = pythag(g,1.0);
				g = d[m] - d[l] + e[l]/(g+sign(r,g));
				s = c = 1.0;
				p = 0.0;
				for (i=m-1;i>=l;i--)
				{
					f = s*e[i];
					b = c*e[i];
					e[i+1] = (r=pythag(f,g));
					if (r == 0.0)
					{
						d[i+1] -= p;
						e[m] = 0.0;
						break;
					}
					s = f/r;
					c = g/r;
					g = d[i+1] - p;
					r = (d[i]-g)*s+2.0*c*b;
					d[i+1] =g+(p=s*r);
					g = c*r-b;
					for (k=0;k<n;k++)
					{
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l)
					continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m != l);
	}
}

vector<EigenReal> mbl::CalcEigen(const matrix<double>& SymmetricMatrix)
{
	const int n = SymmetricMatrix.NrCols();

	double **A = conversion_matrix<double>(SymmetricMatrix);
	double *d = new double [n];
	double *e = new double [n];
	//vector<double> d(n),e(n);

	HouseHolderReduction(A,d,e,n);
	QL_algo(d,e,A,n);

	vector<EigenReal> eigen(n);
	vector<double> u(n);
	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
			u[j] = A[j][i];
		if (accumulate(u.begin(),u.end(),0.0) < 0.0)
			u = -u;
		eigen[i] = EigenReal(d[i],u);
	}
	freematrix<double>(A,n);
	delete[] d;
	delete[] e;

	sort(eigen.begin(),eigen.end(),compare_EigenReal);

	return eigen;
}
