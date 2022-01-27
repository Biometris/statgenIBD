#ifndef BINARYFILE_HEADER
#define BINARYFILE_HEADER

#include <iostream>
#include <fstream>
#include "util_genetics.h"
#include "ibdexcept.h"
#include "misc.h"

// member functions write and read defined for the fundamental types:
// 
// char, signed char, unsigned char
// short, unsigned short
// int, unsigned int
// long, unsigned long
// float 
// double 
// long double
// bool

namespace ibd
{

class oBinFile
{
public:
	oBinFile() {} 
	oBinFile(const std::string& filename) { open(filename);} 
    ~oBinFile() { close(); }
	
	void open(const std::string& filename) 
	{  
		outp.open(filename.c_str(),std::ios::binary);
		if (!outp)
			throw ibd_error("Cannot open file " + filename);
	}
	void close() { outp.close(); }
	
	void write(char * x, size_t size) { outp.write(x,size); }

	// write functions for fundamental types:
	void write(char x)			 { write((char *)&x, sizeof(char)); }
	void write(unsigned char x)  { write((char *)&x, sizeof(unsigned char)); }
	void write(signed char x)    { write((char *)&x, sizeof(signed char)); }
	void write(short x)          { write((char *)&x, sizeof(short)); }
	void write(unsigned short x) { write((char *)&x, sizeof(unsigned short)); }
	void write(int x)			 { write((char *)&x, sizeof(int)); }
	void write(unsigned int x)   { write((char *)&x, sizeof(unsigned int)); }
	void write(long x)			 { write((char *)&x, sizeof(long)); }
	void write(unsigned long x)  { write((char *)&x, sizeof(unsigned long)); }
	void write(float x)			 { write((char *)&x, sizeof(float)); }
	void write(double x)		 { write((char *)&x, sizeof(double)); }
	void write(long double x)	 { write((char *)&x, sizeof(long double)); }
    void write(bool x)			 { write((char *)&x, sizeof(bool)); }
private:
	std::ofstream outp;
};

// save element x of fundamental type
inline oBinFile& operator<<(oBinFile& os, char x)	        { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, signed char x)    { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, unsigned char x)  { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, short x)          { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, unsigned short x) { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, int x)            { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, unsigned int x)   { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, long x)           { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, unsigned long x)  { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, float x)          { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, double x)         { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, long double x)    { os.write(x);return os; }
inline oBinFile& operator<<(oBinFile& os, bool x)           { os.write(x);return os; }


template <class T1, class T2>
oBinFile& operator<<(oBinFile& os, const std::pair<T1,T2>& data)
{
	os << data.first << data.second;
	return os;
}

inline oBinFile& operator<<(oBinFile& os, const std::string& data)
{  
	long int N = data.size();
	os.write(N);
	for (std::string::const_iterator it=data.begin();it != data.end(); it++)
		os << *it;
	return os;
}
 
template <class T>
oBinFile& operator<<(oBinFile& os, const std::vector<T>& data)
{  
	long int N = data.size();
	os.write(N);
	for (typename std::vector<T>::const_iterator it=data.begin();it != data.end(); it++)
		os << *it;
	return os;
}

template <class K, class V>
oBinFile& operator<<(oBinFile& os, const std::map<K,V>& data)
{
	long int N = data.size();
	os.write(N);
	for (typename std::map<K,V>::const_iterator it = data.begin(); it != data.end(); it++)
		os << it->first << it->second;
	return os;
}


class iBinFile
{
public:
	iBinFile() {} 
	iBinFile(const std::string& filename) { open(filename);} 
    ~iBinFile() { close(); }
	
	void open(const std::string& filename) 
	{  
		inp.open(filename.c_str(),std::ios::binary);
		if (!inp)
			throw ibd_error("Cannot open file " + filename);
	}
	void close() { inp.close(); }
	
	void read(char * x, size_t size) { inp.read(x,size); }

	// read functions for fundamental types:
	void read(char& x)			 { read((char *)&x, sizeof(char)); }
	void read(unsigned char& x)  { read((char *)&x, sizeof(unsigned char)); }
	void read(signed char& x)    { read((char *)&x, sizeof(signed char)); }
	void read(short& x)          { read((char *)&x, sizeof(short)); }
	void read(unsigned short& x) { read((char *)&x, sizeof(unsigned short)); }
	void read(int& x)			 { read((char *)&x, sizeof(int)); }
	void read(unsigned int& x)   { read((char *)&x, sizeof(unsigned int)); }
	void read(long& x)			 { read((char *)&x, sizeof(long)); }
	void read(unsigned long& x)  { read((char *)&x, sizeof(unsigned long)); }
	void read(float& x)		     { read((char *)&x, sizeof(float)); }
	void read(double& x)		 { read((char *)&x, sizeof(double)); }
	void read(long double& x)	 { read((char *)&x, sizeof(long double)); }
	void read(bool& x)           { read((char *)&x, sizeof(bool)); }

private:
	std::ifstream inp;
};

// read element x of fundamental type
inline iBinFile& operator>>(iBinFile& is, char& x)	         { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, signed char& x)    { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, unsigned char& x)  { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, short& x)          { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, unsigned short& x) { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, int& x)            { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, unsigned int& x)   { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, long& x)           { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, unsigned long& x)  { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, float& x)          { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, double& x)         { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, long double& x)    { is.read(x);return is; }
inline iBinFile& operator>>(iBinFile& is, bool& x)			 { is.read(x);return is; }

template <class T1, class T2>
iBinFile& operator>>(iBinFile& is, std::pair<T1,T2>& data)
{
	is >> data.first >> data.second;
	return is;
}

inline iBinFile& operator>>(iBinFile& is, std::string& data)
{  
	long int N;
	is.read(N);
	data.resize(N);
	for (std::string::iterator it=data.begin();it != data.end(); it++)
		is >> *it;
	return is;
}

template <class T>
iBinFile& operator>>(iBinFile& is, std::vector<T>& data)
{  
	long int N;
	is.read(N);
	data.resize(N);
	for (typename std::vector<T>::iterator it=data.begin();it != data.end(); it++)
		is >> *it;
	return is;
}

template <class K, class V>
iBinFile& operator>>(iBinFile& is, std::map<K,V>& data)
{
	long int N;
	is.read(N);
	data.clear(); // niet helemaal zeker of dit moet?!
	for (long int i=0;i<N;i++)
	{
		std::pair<K,V> x;
		is >> x.first >> x.second;
		data.insert(x);
	}
	return is;
}


}

#endif
