#include "initialize.hpp"
#include "rarray.hpp"

//Initializes cell array with a specified fraction of alive cells
void initialize_cells(rvector<bool>& alive_status, double target_fraction) {
    double fill = 0.0;

    //Iterate over each cell
    for (bool& alive : alive_status) {
        //Increment fill
        fill += target_fraction;

        //If accumlation gets to 1 or larger set cell as alive, otherwise keep cell dead
        if (fill >= 1.0) {
            alive = true;
            fill -= 1.0;
        } else {
            alive = false;
        }
    }
}