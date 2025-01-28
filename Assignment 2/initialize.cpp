#include "initialize.hpp"
#include "rarray.hpp"

void initialize_cells(rvector<bool>& alive_status, double target_fraction) {
    double fill = 0.0;
    for (bool& alive : alive_status) {
        fill += target_fraction;
        if (fill >= 1.0) {
            alive = true;
            fill -= 1.0;
        } else {
            alive = false;
        }
    }
}