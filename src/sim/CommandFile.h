#ifndef COMMAND_FILE_HEADER
#define COMMAND_FILE_HEADER

#include <map>
#include <string>
#include <set>
#include <iostream>
#include <sstream>
#include "util_genetics.h"
#include "markerscore.h"
#include "misc.h"

namespace ibd
{

typedef 
struct 
{
	int line_nr;
	std::string argument;
	std::string filename;
	void error(const std::string& str) const
		{ throw ibd_error(filename + "(" + itostr(line_nr) + "): " + str);}
} line_command;

typedef std::map<std::string, std::pair<int,int> > DefCommands; // set of commands

class Commands : public std::multimap<std::string,line_command>
{
public:
	Commands() { }
	int GetMaxNr(const std::string& command) const;
	void check(const std::string& command, bool multi_command) const;

private:
	std::map<std::string, int> max_commands; 
	friend void read_command_file(Commands& commands, const DefCommands& command_set,
					   const std::string& filename);
};

inline void AddCommand(DefCommands& def_commands, const std::string& command)
{ def_commands[command] = std::pair<int,int>(1,1); }

inline void AddOptionalCommand(DefCommands& def_commands, const std::string& command)
{ def_commands[command] = std::pair<int,int>(0,1); }

inline void AddMultiCommand(DefCommands& def_commands, const std::string& command, 
							int min, int max)
{ def_commands[command] = std::pair<int,int>(min,max); }

const line_command& GetCommand(const Commands& commands, const std::string& str);

Commands read_command_file(const DefCommands& command_set, const std::string& filename);
void read_command_file(Commands& commands, const DefCommands& command_set,
					   const std::string& filename);

int test_command_file();

template <class T>
bool read(T& val, const Commands& commands, const std::string& command)
{
	commands.check(command,false);
	Commands::const_iterator it = commands.find(command);		
	if (it == commands.end()) return false;
	const line_command& lc = it->second;
	std::istringstream line_stream(lc.argument);
	line_stream >> val;
	if (!line_stream.eof() | line_stream.bad())
		lc.error("wrong arguments");
	return true;
}

template <class T1,class T2>
bool read(T1& val1, T2& val2, const Commands& commands, const std::string& command)
{
	commands.check(command,false);
	Commands::const_iterator it = commands.find(command);
	if (it == commands.end()) return false;

	const line_command& lc = it->second;
	std::istringstream line_stream(lc.argument);
	line_stream >> val1 >> val2;
	if (!line_stream.eof() | line_stream.bad())
		lc.error("wrong arguments");
	return true;
}

template <class T1,class T2, class T3>
bool read(T1& val1, T2& val2, T3& val3, const Commands& commands, 
										const std::string& command)
{
	commands.check(command,false);
	Commands::const_iterator it = commands.find(command);
	if (it == commands.end()) return false;

	const line_command& lc = it->second;
	std::istringstream line_stream(lc.argument);
	line_stream >> val1 >> val2 >> val3;
	if (!line_stream.eof() | line_stream.bad())
		lc.error("wrong arguments");
	return true;
}

template <class T1,class T2, class T3, class T4>
bool read(T1& val1, T2& val2, T3& val3, T4& val4, const Commands& commands, 
										const std::string& command)
{
	commands.check(command,false);
	Commands::const_iterator it = commands.find(command);
	if (it == commands.end()) return false;

	const line_command& lc = it->second;
	std::istringstream line_stream(lc.argument);
	line_stream >> val1 >> val2 >> val3 >> val4;
	if (!line_stream.eof() | line_stream.bad())
		lc.error("wrong arguments");
	return true;
}

template <class T>
bool multi_read(std::vector<T>& x, const Commands& commands, const std::string& command)
{
	commands.check(command,true);
	typedef Commands::const_iterator Iter; 
	std::pair<Iter,Iter> range = commands.equal_range(command);
	if (range.first == range.second)
		return false;
	T value;
	x.clear();
	for (Iter iter = range.first; iter != range.second; ++iter)
	{
		const line_command& lc = iter->second;
		std::istringstream line_stream(lc.argument);
		line_stream >> value;
		x.push_back(value);
		if (!line_stream.eof() || line_stream.bad())
			lc.error("wrong arguments");
	}
	return true;
}

}

#endif
