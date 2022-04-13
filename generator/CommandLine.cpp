/*! \file
 *
 * \brief Some helpers for parsing the command line (source)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */

#include <string>
#include <stdexcept>
#include "generator/CommandLine.hpp"
#include "generator/Types.hpp"


// For options parsing
std::string GetNextArg(int & i, int argc, char ** argv)
{
    if(i >= argc)
        throw std::runtime_error("Error - no more arguments!");

    return argv[i++];
}


int GetIArg(int & i, int argc, char ** argv)
{   
    std::string str = GetNextArg(i, argc, argv);
    try {

        return stoi(str);
    }
    catch(...)
    {
        std::string errstr("Cannot convert to int: ");
        errstr += str;
        throw std::runtime_error(errstr);
    }
}


std::string GetNextArg(size_t & i, const std::vector<std::string> & opt)
{
    if(i >= opt.size())
        throw std::runtime_error("Error - no more arguments!");

    return opt[i++];
}


int GetIArg(size_t & i, const std::vector<std::string> & opt)
{
    std::string str = GetNextArg(i, opt);
    try {
        return stoi(str);
    }
    catch(...)
    {
        std::string errstr("Cannot convert to int: ");
        errstr += str;
        throw std::runtime_error(errstr);
    }
}


std::vector<std::string> ParseCommonOptions(OptionMap & options, int argc, char ** argv)
{
    std::vector<std::string> ret;

    int i = 1;
    while(i < argc)
    {
        std::string argstr(GetNextArg(i, argc, argv));
        if(argstr == "-ve")
            options[Option::ExternalVRR] = GetIArg(i, argc, argv);
        else if(argstr == "-vg")
            options[Option::GeneralVRR] = GetIArg(i, argc, argv);
        else if(argstr == "-he")
            options[Option::ExternalHRR] = GetIArg(i, argc, argv);
        else if(argstr == "-hg")
            options[Option::GeneralHRR] = GetIArg(i, argc, argv);
        else if(argstr == "-p")
            options[Option::FinalPermute] = 1;
        else if(argstr == "-omp")
            options[Option::OpenMPTarget] = 1;
        else
            ret.push_back(argstr);
    }

    return ret;
}


