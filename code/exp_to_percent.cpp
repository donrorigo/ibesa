#include <iostream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;

int main(int argc, char const *argv[])
{
    

    ifstream proteins;
    string amino_sequence;
    double hypervolume;
    string code;
    

    proteins.open(argv[1]);
    while(!proteins.eof())
    {
        getline(proteins, code, ' ');
        hypervolume = stod(code);
        std::cout << abs(hypervolume)*100 << std::endl;
    }
    proteins.close();

    return 0;
}
