#include "output.hpp"
#include "rarray.hpp"

void output(const rvector<bool>& alive_status, int step, char alive_char, char dead_char) {
    // Calculate the fraction of alive cells
    double fraction = std::count(alive_status.begin(), alive_status.end(), true) / double(alive_status.size());

    // Create the string representation of the cells
    std::string string_representation(alive_status.size(), ' ');
    for (int i = 0; i < alive_status.size(); i++) {
        if (alive_status[i]) {
            string_representation[i] = alive_char;
        } else {
            string_representation[i] = dead_char;
        }
    }

    // Print the output
    std::cout << step << "\t" << string_representation << " " << fraction << "\n";
}

