#include <fstream>
#include <sstream>
#include "util_genetics.h"
#include "CommandFile.h"

using namespace std;
using namespace ibd;

int Commands::GetMaxNr(const string& command) const
{
    map<string,int>::const_iterator p = max_commands.find(command);
    if (p == max_commands.end()) return 0;
    return p->second;
}

void Commands::check(const string& command, bool multi_command) const
{
	const int M = GetMaxNr(command);
	if (M == 0)
		throw ibd_error("Command " + command + " not defined");
	if (multi_command)
	{
		if (M <= 1)
			throw ibd_error("Use function read() for command " + command);
	}
	else
	{
		if (M > 1)
			throw ibd_error("Use function multi_read() for command " + command);
	}
}


ibd::Commands ibd::read_command_file(const DefCommands& command_set, const string& filename)
{
	Commands commands;
	read_command_file(commands,command_set,filename);
	return commands;
}

void ibd::read_command_file(Commands& commands, const DefCommands& command_set,
							const string& filename)
{
	ifstream inp;
	OpenFile(inp,filename);

	int line_number = 0;
	string line;

	line_command lc;
	lc.filename = filename;

	while (getline(inp,line))
	{
		line_number++;
		line = remove_comment(line);
		if (!line.empty())
		{
			string command;
			lc.line_nr = line_number;
			istringstream line_stream(line);
			line_stream >> command;
			map<string, pair<int,int> >::const_iterator it = command_set.find(command);
			if (it == command_set.end())
				lc.error("unknown command " + command);
			getline(line_stream,lc.argument);
			commands.insert(pair<string,line_command>(command,lc));
			if (it->second == pair<int,int>(1,1) && commands.count(command) > 1)
				lc.error("multiple definition of command " + command);
		}
	}
	typedef map<string, pair<int,int> >::const_iterator Iter;
	for (Iter it  = command_set.begin(); it != command_set.end(); ++it)
	{
		const string& command = (*it).first;
		const int& min = (*it).second.first;
		const int& max = (*it).second.second;
		const int& cnt = commands.count(command);
		if (cnt == 0 && min > 0)
			throw ibd_error("Command " + command + " not defined");
		if (cnt < min || cnt > max)
		{
			string error = "Wrong number of command " + command + " in file " + filename;
			throw ibd_error(error);
		}
		commands.max_commands[command] = max;
	}
}

const ibd::line_command& ibd::GetCommand(const Commands& commands, const string& str)
{
	Commands::const_iterator it = commands.find(str);
	if (it != commands.end())
		return it->second;
	else
		throw ibd_error("function GetCommand(): Undefined command " + str);
}


/* Test of script files (file with commands)
This function uses the file "script.txt" as an example:

CHILD	  Piet
PHONENUMB 426288 477316
CHILD	  Kees
CHILD	  Clara
LASTNAME  Klaasen
FIRSTNAME Jan
POSTCODE  6700 AB
BIRTHDATE 5 April 1970

*/
// int ibd::test_command_file()
// {
// 	// Definition of the commands in the script file. For the script file we
// 	// require one unique definition of LASTNAME, POSTCODE, BIRTHDATE, PHONENUMB,
// 	// the command FIRSTNAME is optional. The minimum number of CHILD commands is
// 	// 0 and the maximum number is 10.
// 	DefCommands defined_commands;
// 	AddCommand(defined_commands,"LASTNAME");
// 	AddCommand(defined_commands,"POSTCODE");
// 	AddCommand(defined_commands,"BIRTHDATE");
// 	AddCommand(defined_commands,"PHONENUMB");
// 	AddOptionalCommand(defined_commands,"FIRSTNAME");
// 	AddMultiCommand(defined_commands,"CHILD",0,10);
//
// 	// read the commands from file script.txt and check for errors
// 	Commands commands = read_command_file(defined_commands,"script.txt");
//
// 	int digits_pc, day, year;
// 	string firstname = "Not defined";
// 	string lastname, char_pc, month;
// 	vector<string> children;
// 	vector<int> phone_numb(2);
//
// 	// read the data in the variables and check for errors
// 	read(firstname,commands,"FIRSTNAME");		 // 1 argument (type string)
// 	read(lastname,commands,"LASTNAME");	         // 1 argument (type string)
// 	read(digits_pc,char_pc,commands,"POSTCODE"); // 2 arguments (types int and string)
// 	read(day,month,year,commands,"BIRTHDATE");   // 3 arguments (types int,string,int)
// 	read(phone_numb,commands,"PHONENUMB");       // 1 argument (vector<int>, size = 2)
//
// 	// this command will read all the CHILD commands
// 	multi_read(children,commands,"CHILD");       // 1 argument (vector of type string)
//
// 	cout << "First name:    " << '\t' << firstname << endl;
// 	cout << "Last name:     " << '\t' << lastname << endl;
// 	cout << "Postcode:      " << '\t' << digits_pc << " " << char_pc << endl;
// 	cout << "Date of birth: " << '\t' << day << " " << month << " " << year << endl;
// 	cout << "Phone numbers: " << '\t' << phone_numb[0] << " " << phone_numb[1] << endl;
// 	cout << "Children:   " << '\t';
// 	if (children.empty())
// 		cout << "None" << endl;
// 	else
// 		for (unsigned int i=0;i<children.size();i++)
// 			cout << children[i] << " ";
// 	cout << endl << endl;
//
// 	return 0;
// }

