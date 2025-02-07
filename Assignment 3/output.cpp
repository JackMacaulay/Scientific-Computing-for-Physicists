// output.cpp
#include "output.h"
#include <fstream>
void writeBinary(const std::string& name,
                 const rvector<double>& a) {
    std::ofstream g(name, std::ios::binary);
    g.write((char*)(&a[0]), a.size()*sizeof(a[0]));
    g.close();
}
void writeText(const std::string& name,
               const rvector<double>& a) {
    std::ofstream f(name);
    for (int i = 0; i < a.size(); i++)
        f << a[i] << std::endl;
    f.close();
}
