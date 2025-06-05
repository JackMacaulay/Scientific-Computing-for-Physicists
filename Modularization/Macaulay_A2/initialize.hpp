#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "rarray.hpp"

//Initializes the array of cells with a target fraction of alive cells.

void initialize_cells(rvector<bool>& alive_status, double target_fraction);

#endif 
