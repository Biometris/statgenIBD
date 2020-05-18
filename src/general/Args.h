/*!
\file
\brief  Reading command-line arguments and options.
\author Martin Boer, Biometris
\date   1998-2007
*/
#ifndef ARGS_HEADER
#define ARGS_HEADER

#include <set>
#include <map>
#include <string>
#include <vector>
#include <iostream>

namespace mbl
{

class Options;
class Args;

class Option
{
public:
	Option(const std::string& k, const std::string& t, const std::string& d);
	std::string PrintType() const;
private:
	std::string key, descr;
	enum Type {NO_ARG, INT_ARG, DOUBLE_ARG, STRING_ARG};
	Type type;
	bool state;

	friend class Options;
	friend class Args;
};


class Options : public std::vector<Option>
{
public:
    Options() {}
    void Add(const std::string& opt, const std::string& descr);
	void print_descriptions(std::ostream& f) const;
	const_iterator find(const std::string& key) const;
};

std::ostream& operator<<(std::ostream& outp, const Options& options);

class Args
{
public:
	Args(const Options& options,
		 std::string args,
		 std::string opt_args,
		 int argc, const char *argv[]);

	Args(std::string format,
		 std::string args,
		 std::string opt_args,
		 int argc, const char *argv[]);

	void init(std::string args, std::string opt_args, int argc, const char *argv[]);

	std::string GetArgument(int i) const { return m_ArgValues[i]; }
	bool GetArgument(const std::string& arg_name, std::string& arg_val) const;
	int GetNrArguments() const { return (int) m_ArgValues.size(); }
	bool GetOption(const std::string& option) const;
	bool GetOption(const std::string& option, int&) const;
	bool GetOption(const std::string& option, double&) const;
	bool GetOption(const std::string& option, std::string&) const;
	bool State() const { return m_state; }
	std::string Error_string() const { return m_error_string; }
	int PrintOptions() const;

private:

	bool error(const std::string& str) { m_error_string = str; return false; }
	bool init_option_format(const std::string& format);
	bool read_args(const std::string& format);
	bool init_arguments(int argc, const char *argv[]);

	Options opts;

	std::set<std::string>			      m_OptVoid;    // init. by *argv[]

	std::map<std::string,int>			  m_OptInt;		// init. by *argv[]
	std::map<std::string,double>		  m_OptDouble;  // init. by *argv[]
	std::map<std::string,std::string>     m_OptString;  // init. by *argv[]

	std::map<std::string,int>			  m_ArgNames;   // init. by format string
	std::vector<std::string>			  m_ArgValues;  // init. by *argv[]
	int	m_min_narg,	m_max_narg;

	// for error handling of constructor
	bool m_state;
	std::string m_error_string;

};

}

#endif
