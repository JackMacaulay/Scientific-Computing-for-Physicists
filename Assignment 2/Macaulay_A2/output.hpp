#ifndef OUTPUT_H
#define OUTPUT_H

#include "rarray.hpp"

//Initializes the array of cells with a target fraction of alive cells.

void output(const rvector<bool>& alive_status, int step, char alive_char, char dead_char);


#endif 
