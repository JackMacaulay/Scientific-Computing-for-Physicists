// analyze.cpp
#include "friction.h"
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

void analyze_data(const std::string& filename,
                  const std::string& outfilename)
{
    std::ifstream fin;
    fin.open(filename);
    long long n=0;
    std::string line;
    while (std::getline(fin,line)) {
        if (line[0] != '#') {
            n++;
        }
    }
    fin.close();
    fin.open(filename);
    rvector<double> t(n);
    rvector<double> z(n);
    n = 0;
    while (std::getline(fin,line) and n < t.size()) {
        if (line[0] != '#') {
            std::stringstream s(line);
            s >> t[n] >> z[n];
            n++;
        }
    }
    fin.close();
    // compute alpha
    double dt = t[1]-t[0];
    rvector<double> dzdt = numdiff(dt, z);
    double alpha = frictionrate(dt, dzdt);
    // report
    std::ostream* fout;
    if (outfilename == "-")
        fout = &std::cout;
    else
        fout = new std::ofstream(outfilename);
    *fout << "# filename " << filename << "\n";
    *fout << "# n " << n << "\n";
    *fout << "# dt " << dt << "\n";
    *fout << "# alpha " << alpha << '\n';
    if (outfilename != "-") {
        delete fout;
        std::cout << "Output written to '" << outfilename << "'.\n";        
    }
}

int read_command_line(int argc, char* argv[],
                      std::string& filename,
                      std::string& outfilename)
{
    using boost::program_options::value;
    boost::program_options::options_description desc("Options for analyze");
    desc.add_options()
        ("help,h",                                      "Print help message")
        ("file,f",    value<std::string>(&filename),    "file from which read the data")
        ("output,o",  value<std::string>(&outfilename), "file to write result (- means console)");
    boost::program_options::variables_map args;
    try {
        store(parse_command_line(argc, argv, desc), args);
        notify(args);
    }
    catch (...) {
        std::cerr << "ERROR in command line arguments!\n" << desc;
        return 2;
    }
    if (args.count("help")) {
        std::cout << "Usage:\n    " << argv[0] << " [OPTIONS]\n" << desc;
        return 1;
    }
    return 0;
}

int main(int argc, char* argv[])
{
    // Define defaults settings
    std::string filename = "testmodel.dat";
    std::string outfilename = "-";
    // Parse the command line
    int status = read_command_line(argc, argv, filename, outfilename);
    if (status > 0) 
        return status - 1;
    // Go to work
    analyze_data(filename, outfilename);
}
