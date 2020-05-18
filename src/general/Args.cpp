/***********************************************************
 *														   *
 *   Martin Boer											   *
 *   Biometris Wageningen-UR & NWO, grant no. 925-01-001    *
 *   July 12, 2002										   *
 *														   *
 *   update october 6, 2006:                                *
 *	save options cout formatting with                      *
 *   ios_base::fmtflags old_options = cout.flags();         *
 *   and restore with cout.flags(old_options);              *
 *                                                          *
 *   see: void Options::print_descriptions procedure        *
 *                                                          *
 *   Dec 24, 2006: Revision of Args and Option classes      *
 ************************************************************/
#include <iomanip>
#include <iostream>
#include <sstream>
#include <Rcpp.h>

#include "Args.h"
#include "convert.h"
#include "mblexcept.h"

using namespace std;
using namespace mbl;

Option::Option(const std::string& k,
               const std::string& t,
               const std::string& d)
  : key(k), descr(d)
{
  state = true;
  if (t == "")
    type = NO_ARG;
  else if (t == "<int>")
    type = INT_ARG;
  else if (t == "<double>")
    type = DOUBLE_ARG;
  else if (t == "<string>")
    type = STRING_ARG;
  else
    state = false;
}

string Option::PrintType() const
{
  if (!state) return "error!";
  switch (type)
  {
  case Option::NO_ARG:	 return "";
  case Option::INT_ARG:	 return "<int>";
  case Option::DOUBLE_ARG: return "<double>";
  case Option::STRING_ARG: return "<string>";
  }
  return "";
}


void Options::Add(const std::string& opt,
                  const std::string& descr)
{
  istringstream i(opt);
  string key,type;
  i >> key >> type;
  push_back(Option(key,type,descr));
}

Options::const_iterator Options::find(const std::string& key) const
{
  const_iterator it;
  for (it=begin();it!=end();it++)
    if (key == it->key) return it;
    return it;
}

void Options::print_descriptions(std::ostream& f) const
{
  // ios_base::fmtflags old_options = cout.flags();
  // f.setf(ios::left);
  for (const_iterator it = begin(); it != end(); ++it)
  {
    const string opt = "-" + it->key;
    const string first = (it == begin()) ? "   options: " : "            ";
    f << first << setw(10) << opt << setw(10) << it->PrintType();
    if (!it->descr.empty())
      f << ": " << it->descr;
    f << endl;
  }
  // cout.flags(old_options);
}

ostream& mbl::operator<<(std::ostream& outp,
                         const Options& options)
{
  options.print_descriptions(outp);
  return outp;
}


void Args::init(std::string args,
                std::string opt_arg,
                int argc,
                const char *argv[])
{
  read_args(args);
  m_min_narg = m_ArgNames.size();
  read_args(opt_arg);
  m_max_narg = m_ArgNames.size();
  m_state = init_arguments(argc,argv);
}

Args::Args(const Options& options,
           std::string args,
           std::string opt_arg,
           int argc,
           const char *argv[])
{
  opts = options;
  init(args,opt_arg,argc,argv);
}

Args::Args(std::string format,
           std::string args,
           std::string opt_arg,
           int argc,
           const char *argv[])
{
  m_state = init_option_format(format);
  if (m_state)
    init(args,opt_arg,argc,argv);
}

int Args::PrintOptions() const
{
  if (!m_state)
    Rcpp::Rcout << "Error: " << m_error_string << endl;
  Rcpp::Rcout << endl << opts << endl << endl;
  return 1;
}

bool Args::GetOption(const std::string& option) const
{
  if (opts.find(option) == opts.end())
    throw mblib_error("Option -" + option + " not defined");

  set<string>::const_iterator it = m_OptVoid.find(option);
  if (it == m_OptVoid.end()) return false;
  return true;
}

template<class T>
bool get_option(const std::string& option,
                const std::map<std::string,T>& map_type,
                const Options& def_options, T& val)
{
  if (def_options.find(option) == def_options.end())
    throw mblib_error("Option -" + option + " not defined");
  typename map<string,T>::const_iterator it = map_type.find(option);
  if (it == map_type.end()) return false;
  val = it->second;
  return true;
}

bool Args::GetOption(const std::string& option, int& value) const
{ return get_option(option,m_OptInt,opts,value); }

bool Args::GetOption(const std::string& option, double& value) const
{ return get_option(option,m_OptDouble,opts,value); }

bool Args::GetOption(const std::string& option, std::string& value) const
{ return get_option(option,m_OptString,opts,value); }

bool Args::GetArgument(const std::string& arg_name, std::string& arg_val) const
{
  map<string,int>::const_iterator it = m_ArgNames.find(arg_name);
  if (it == m_ArgNames.end())
    throw mblib_error("GetArgument: argument '" + arg_name + "' not defined");
  if (it->second >= (int) m_ArgValues.size())
    return false;
  arg_val = m_ArgValues[it->second];
  return true;
}

bool Args::read_args(const std::string& format)
{
  std::istringstream form(format);
  string name;
  while (form >> name)
  {
    if (m_ArgNames.find(name) != m_ArgNames.end())
      return error("redefinition of argument " + name);
    int size =  m_ArgNames.size();
    m_ArgNames[name] = size;
  }
  return true;
}

bool Args::init_option_format(const std::string& format)
{
  char c;
  string arg_type;
  std::istringstream form(format);
  while (form >> c)
  {
    if (c != '-')
      return error("Error while reading format string for options");
    string option_name,option_type;
    form >> option_name;
    if (opts.find(option_name) != opts.end())
      return error("Multiple occurence of option -" + option_name + " in format");
    form >> c;
    form.putback(c);
    if (c=='-' || c==EOF)
      arg_type = "";
    else
    {
      form >> option_type;
      if (option_type == "<int>")
        arg_type = option_type;
      else if (option_type == "<double>")
        arg_type = option_type;
      else if (option_type == "<string>")
        arg_type = option_type;
      else
      {
        string str = "Unknown argument " + option_type
        + " for option " + option_name;
        return error(str);
      }
    }
    opts.Add(option_name + " " + option_type,"");
  }
  return true;
}


bool Args::init_arguments(int argc, const char *argv[])
{
  for (Options::const_iterator it=opts.begin();it!=opts.end();it++)
    if (!it->state)
      return error("Error in definition options");

    int nargs=0;
    set<string> SetOptions;
    for (int i=1;i<argc;i++)
    {
      istringstream inp(argv[i]);
      char c;
      inp >> c;
      if (c == '-') // option
      {
        string option;
        inp >> option;
        if (SetOptions.find(option) != SetOptions.end())
          return error("Multiple occurence of option -" + option);
        else
          SetOptions.insert(option);

        Options::const_iterator it = opts.find(option);
        if (it == opts.end())
          return error("Unknown option -" + option);
        if (it->type == Option::NO_ARG)
          m_OptVoid.insert(option);
        else
        {
          i++;
          if (i == argc)
            return error("No argument for option -" + option);
          string val(argv[i]);
          switch (it->type)
          {
          case Option::INT_ARG:
            m_OptInt[option] = convertTo<int>(val); break;
          case Option::DOUBLE_ARG:
            m_OptDouble[option] = convertTo<double>(val); break;
          case Option::STRING_ARG:
            m_OptString[option] = val; break;
          case Option::NO_ARG:
            break;
          }
        }
      }
      else
      {
        nargs++;
        if (nargs > m_max_narg)
          return error("Too many arguments");
        m_ArgValues.push_back(argv[i]);
      }
    }
    if (nargs < m_min_narg)
      return error("Too few arguments");
    return true;
}

