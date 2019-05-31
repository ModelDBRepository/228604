#!/usr/bin/python

import imp
import sys
import string
from os.path import dirname
import json

headerstr = """
#ifndef AUTOPARAMS_H_
#define AUTOPARAMS_H_

// automatically generated. change 'parameters' instead

#include <iostream>
#include <string>

using namespace std; // for string

namespace params 
{

$param_declarations

	int acquire(int argc, char** argv);
	void log(std::ostream& out);
};

#endif // AUTOPARAMS_H_
"""

sourcestr = """
// automatically generated. change 'parameters' instead

#include <boost/program_options.hpp>
#include <iostream>
#include "autoparams.h"

namespace po = boost::program_options;

namespace params
{

$param_definitions

}


int params::acquire(int argc, char** argv)
{
	po::options_description od("options");
	po::variables_map options;
	od.add_options()
		("help", "show this help")

$options

		;

	po::store(po::command_line_parser(argc, argv) 
		.options(od)
		.style(po::command_line_style::default_style ^ po::command_line_style::allow_guessing)
		.run(), options);
	po::notify(options);
	if (options.count("help"))
        {
                std::cout << od;
		return 1;
        }
	else return 0;
}

void params::log(std::ostream& out)
{

$logging	

}
"""

params_def_filename = sys.argv[1]

header_out_fn = "./autoparams.h"
source_out_fn = "./autoparams.cpp"

def c_decl(p_entry):
    name, ctype, defval = p_entry
    return "extern " + ctype + " " + name + ";\n"

def c_def(p_entry):
    name, ctype, defval = p_entry
    if ctype == "string":
        return ctype + " " + name + '("' + str(defval) + '");\n'
    else:
        return ctype + " " + name + " = " + str(defval) + ";\n"

try:
    fh = open(params_def_filename, "r")
    prms = json.load(fh)
    fh.close()
except IOError:
    print("ERROR: Cannot load " + params_def_filename)
    sys.exit(1)

declstr = ""
defstr = ""
optionstr = ""
loggingstr = ""

for k in sorted([k for k in prms.keys() if not k.startswith("_")]):
    if not (hasattr(prms[k], "has_key") and prms[k].has_key("c_type")):
        # try inferring the type. works only for int, double, string
        val = prms[k]
        t = type(val)
        tmap = {type(""): "string",
                type(u""): "string",
                type(1): "double", # better safe than sorry
                type(1.5): "double"}
        try:
            ctype = tmap[t]
        except KeyError:
            print("ERROR: cannot infer type of '%s'" % k)
            print("\t('%s' is not in %s)" % (t, str(tmap.keys())))
            sys.exit(1)
    else:
        ctype = prms[k]["c_type"]
        val = prms[k]["value"]
    p = (k, ctype, val)
    name, ctype, defval = p
    declstr += c_decl(p)
    defstr += c_def(p)
    optionstr += '("' + name + '", po::value<' + str(ctype) + '>(&params::' +  name + '))\n'
    loggingstr += 'out << "' + name + ' = " << params::' + name + ' << std::endl;\n'
        
ht = string.Template(headerstr)
st = string.Template(sourcestr)

fh = open(header_out_fn, "w")
fh.write(ht.substitute({"param_declarations": declstr}))
fh.close()

fh = open(source_out_fn, "w")
fh.write(st.substitute({"param_definitions": defstr, "options": optionstr, "logging": loggingstr}))
fh.close()



