#ifndef TIME_STEP_H
#define TIME_STEP_H

#include "rarray.hpp"

//Initializes the array of cells with a target fraction of alive cells.

void time_step(rvector<bool>& alive_status);

#endif 
