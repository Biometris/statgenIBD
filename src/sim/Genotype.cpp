#include <iostream>
#include <iomanip>
#include "Genotype.h"

using namespace std;

ostream& ibd::operator<<(ostream& outp, const Genotype& g)
{
	int str_length = 2;
	int nr_whitespaces = outp.width() - str_length;
	outp.width(0);
	for (int i=0;i<nr_whitespaces;i++)
		outp << ' ';
	g.print(outp);
	return outp;
}

int ibd::compare(const Genotype& g1, const Genotype& g2)
{
	if (g1.x < g2.x) return -1;
	if (g1.x > g2.x) return  1;

	if (g1.y < g2.y) return -1;
	if (g1.y > g2.y) return  1;
	
	return 0;
}	
