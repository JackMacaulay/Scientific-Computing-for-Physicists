#include <iostream>
#include <cmath>
#include <rarray>
#include <string>
#include <cstdlib>
#include <random>

//Structure to hold simulation parameters
struct SimulationParams{
	int M;
	int N;
	double p;
	double g;
	int K;
	int S;
};

//Parsing simulation parameters
SimulationParams parse_arguments(int argc, char* argv[]){
	if (argc != 7){
		std::cerr <<"Specify all parameters: " << argv[0] << " M, N, p, g, K, S \n";
		std::exit(1);
	}
	SimulationParams params;
	params.M = std::stoi(argv[1]);
	params.N = std::stoi(argv[2]);
	params.p = std::stod(argv[3]);
	params.g = std::stod(argv[4]);
	params.K = std::stoi(argv[5]);
	params.S = std::stoi(argv[6]);

	return params;
}


//Generates a 2D lattice of dimensions N x M
rarray<int, 2> generateLattice(const SimulationParams &params, std::mt19937 &gen) {
	//Initialize lattice
	rarray<int, 2> lattice(params.N, params.M);

	//Create a Bernoulli distribution used to populate cells
	std::bernoulli_distribution dist(params.p);

	//Populating lattice cells
	for (int i=0; i < params.N; ++i) {
		for (int j=0; j < params.M; ++j) {
			lattice[i][j] = dist(gen) ? 1 : 0;
		}
	}
	return lattice;
}


//Simulates one walker starting from an empty cell in the top row
bool simulateWalker(const rarray<int, 2>& lattice, int start_col, const SimulationParams &params, std::mt19937 &gen) {

	int M = params.M, N = params.N;
	int maxSteps =  params.S * (M * M + N * N);
	int row =  0, col = start_col;
	int steps = 0;

	//Main simulation loop - ends once limit S(M^2 + N^2) is reached or walker reaches bottom
	while (steps < maxSteps) {
		//Check if walker has reached the bottom
		if (row == M - 1)
			return true;

		std::vector<double> weights;
		std::vector<std::pair<int, int>> moves;
		//Check for valid left move
		if (col > 0 && lattice[row][col - 1] == 1) {
			moves.push_back({0, -1});
			weights.push_back(1.0);
		}
		//Check right move
		if (col < M - 1 && lattice[row][col+1] == 1) {
			moves.push_back({0, 1});
			weights.push_back(1.0);
		}
		//Check upward move
		if (row > 0 && lattice[row - 1][col] == 1) {
			moves.push_back({-1, 0});
			weights.push_back(1.0 / params.g);
		}
		//Check downward move
		if (row < N - 1 && lattice[row + 1][col] == 1) {
			moves.push_back({1, 0});
			weights.push_back(params.g);
		}
		//If no moves available walker is trapped
		if (moves.empty()) {
			return false;
		}

		//Using a discrete distribution to choose a move based on the weights
		std::discrete_distribution<> dist(weights.begin(), weights.end());
		int index = dist(gen);

		//Update the walker's position
		row += moves[index].first;
		col += moves[index].second;

		++steps;
	}
	return false;
}



int main(int argc, char* argv[]) {
	SimulationParams params = parse_arguments(argc, argv);

	//Initialize random number generator with fixed seed
	std::mt19937 gen(12345);

	rarray<int, 2> lattice = generateLattice(params, gen);

	int totalWalkers = 0;
	int successfulWalkers  = 0;

	//Main loops
	// For each column in top row launch walkers if cell is empty
	for (int col = 0; col < params.M; ++col) {
		if (lattice[0][col] == 1) {
			//Launch walkers and if it reaches bottom count as success
			for (int i = 0; i < params.K; ++i) {
				++totalWalkers;
				if (simulateWalker(lattice, col, params, gen))
					++successfulWalkers;
			}
		}
	}

	//Calculate number of succesful walkers and print results
	double fraction = (totalWalkers > 0) ? static_cast<double>(successfulWalkers) / totalWalkers: 0;

	std::cout << "Total Walkers: " << totalWalkers << "\n"
		  << "Successful walkers: " << successfulWalkers << "\n"
		  << "Fraction of successful walkers: " << fraction << "\n";
	return 0;
}

