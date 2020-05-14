/*!
\file 
\brief  Definition of classes ParentsInd and Pedigree.
\author Martin Boer, Biometris
\date   2006-2010
*/
#ifndef PEDIGREE_PARENTSIND_HEADER
#define PEDIGREE_PARENTSIND_HEADER

#include <set>
#include <vector>
#include <string>

#include "IndProp.h"

//! Defines the parents of an individual 
class ParentsInd : private std::pair<int,int>
{
public:
	//! Default constructor.
	/*! 
	    Undefined parents, i.e. we assume that the individual is an inbred founder
	*/
	ParentsInd() : std::pair<int,int>(-1,-1) {} 

	//! Constructor, defines the parents of an individual.
	/*! 
	  \param par1 first parent
      \param par2 second parent 
	*/
	ParentsInd(int par1, int par2) : std::pair<int,int>(par1,par2) {} 

	//! Returns true if the individual is a founder
	bool IsFounder() const { return this->first < 0; } 

	//! Returns first (x=0) or second (x=1) parent of an individual
	const int& operator[](int x) const;
};

//! vector with parents of the ordered individuals
/*
   Class Pedigree is used in the calculations of IBD-probabilities. 
   It consists of an ordered list of individuals. 
*/
class Pedigree : public std::vector<ParentsInd>
{
public:
	//! Default constructor: defines an empty pedigree
	Pedigree() {}

	//! Constructor
	Pedigree(const std::vector<IndProp>& pop);

	//! Constructor
	Pedigree(const std::vector<ParentsInd>& vecparents);
	
	//! Number of founders in pedigree
	int Nfounders() const { return Nfnd; }

	//! Number of non-founders in pedigree
	int Nnonfounders() const { return this->size() - Nfnd; }
	
	//! Returns a vector of indices of all the founders in pedigree
	std::vector<int> GetNdxFnd() const; 
private:
	int Nfnd;  /*!< Number of founders */
};

#endif


