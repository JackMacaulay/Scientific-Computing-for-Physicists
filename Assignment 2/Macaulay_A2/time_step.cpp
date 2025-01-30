#include "rarray.hpp"
#include "time_step.hpp"

//Function to enforce cell survival and birth rules
bool is_cell_alive_next(const rvector<bool>& cell_alive, int index)
{
    int num_cells = cell_alive.size();
    int alive_neighbours = 0;

    //Iterate through neighbours within 2 cells
    for (int offset = -2; offset <= 2; offset++)
        if (offset != 0 and cell_alive[(index+offset+num_cells)%num_cells]) // modulus (%) enforces periodic boundary conditions
            alive_neighbours++;

    //Applying survival and birth rules
    if (cell_alive[index] and
        (alive_neighbours == 2 or alive_neighbours == 4))
        return true;
    else if (not cell_alive[index] and
             (alive_neighbours == 2 or alive_neighbours == 3))
        return true;
    else
        return false;
}

void time_step(rvector<bool>& alive_status) {
    rvector<bool> next_alive_status(alive_status.size());

    // Update the state of each cell
    for (int i = 0; i < alive_status.size(); i++) {
        next_alive_status[i] = is_cell_alive_next(alive_status, i);
    }

    // Swap current state with the new state
    std::swap(alive_status, next_alive_status);
}