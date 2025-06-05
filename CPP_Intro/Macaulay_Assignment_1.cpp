#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
using namespace std;


bool readDataFromFile(const string& filename, vector<double>& time, vector<double>& x, vector<double>& y);
bool saveResultsToFile(const string& filename, const vector<double>& timeRollAvg, const vector<double>& normRollAvg);
vector<double> computeNorm(const vector<double>& x, const vector<double>& y);
vector<double> computeRollingAverage(const vector<double>& data, int windowSize);


int main(){

    //Opening data file relevant to assignment
    ifstream inputFile("in.dat");

    //Vectors to store the three columns of data
    vector<double> time;
    vector<double> x;
    vector<double> y;

    //Getting filename from user
    string filename;
    string outputfile;
    int n;
    cout << "Enter the name of the data file: ";
    cin >> filename;
    cout << "Specify the width of the running average: ";
    cin >> n;
    cout << "Specify the name of the output file: ";
    cin >> outputfile;

    //Opening the file and populating aforementioned vectors
    if (readDataFromFile(filename, time, x, y)){
        
        //Calculating norm and rolling averages
        vector<double> norms = computeNorm(x,y);
        vector<double> timeRollAvg = computeRollingAverage(time, n);
        vector<double> normRollAvg = computeRollingAverage(norms, n);
        
        // Save results to the specified file
        if (saveResultsToFile(outputfile, timeRollAvg, normRollAvg)) {
            cout << "Results successfully saved to '" << outputfile << "'!" << endl;
        } else {
            cout << "Failed to save results to file." << endl;
        }
    } else {
        cout <<"File ain't there...";
    }

    return 0;
}


// Reads data from a file into three separate vectors: time, x, and y
// Returns true if reading is successful, false otherwise
bool readDataFromFile(const string& filename, vector<double>& time, vector<double>& x, vector<double>& y){
    ifstream infile(filename);
    if (!infile){
        return false;
    }

    double v1, v2, v3;
    // Read file line-by-line, extracting three values at a time
    while(infile >> v1 >> v2 >> v3){
        time.push_back(v1);
        x.push_back(v2);
        y.push_back(v3);
    }

    infile.close();
    return true;
}


// Computes the norm for corresponding elements in x and y
vector<double> computeNorm(const vector<double>& x, const vector<double>& y){
    vector<double> norms;

    // Calculate the norm for each (x, y) pair
    for (int i = 0; i < x.size(); ++i) {
        norms.push_back(sqrt(x[i] * x[i] + y[i] * y[i]));
    }
    return norms;
}

// Computes the rolling average of a data vector given a window size
vector<double> computeRollingAverage (const vector<double>& data, int windowSize){
    vector<double>rollAvg;
    
    int n = data.size();

    // Compute the average over the window size
    for (int i=0; i < n; ++i){
        
        // Calculating start and end indicies
        int start = max(0, i - (windowSize - 1));
        int end   = i;

        // Calculate the sum of the elements in the window
        double sum = 0;
        int count = 0;
        for (int j = start; j <= end; ++j){
            sum += data[j];
            count++;
        }

        // Calculate average
        rollAvg.push_back(sum / count);
    }

    return rollAvg;
    
}

//Saves rolling average results to file
bool saveResultsToFile(const string& filename, const vector<double>& timeRollAvg, const vector<double>& normRollAvg) {
    ofstream outfile(filename);
    if (!outfile) {
        return false;
    }
      
    outfile << "# Time_Avg Norm_Avg\n";

    // Write data
    for (size_t i = 0; i < timeRollAvg.size(); ++i) {
        outfile << timeRollAvg[i] << " " << normRollAvg[i] << "\n";
    }

    outfile.close();
    return true;
}
