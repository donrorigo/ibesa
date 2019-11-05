#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <vector>
#include <string>
#include <iomanip>
#include <boost/tokenizer.hpp>

using namespace std;


int main(int argc, char const *argv[])
{
    if(argc != 4) 
    {
        cout << "El formato debe ser: ./" << argv[0] << " fichero.txt 'codigo de la proteina' 'numero de normalizado' " << endl;
        return 1;
    }
    int j;
    double size;
    string file = argv[1], sequence, code = argv[2];
    float number = std::stof(argv[3]);

    std::vector<std::string> result;
    ofstream output("NORMALIZADO: " + code + "_" +std::to_string(number));
    ifstream input;
    output << fixed;
    output << setprecision(15);
    input.open(file);
    while(!input.eof())
    {
        getline(input, sequence);
        boost::tokenizer<boost::char_separator<char>> tokens(sequence, boost::char_separator<char>(" "));
        std::vector<std::string> result(tokens.begin(), tokens.end());
        for(int i=0; i<result.size(); i++){
            switch (i)
            {
            case 0:
                output << result[i] << " ";
                break;
            case 1:
                output << (stod(result[i])/number) << " ";
                break;
            case 2:
                output << result[i] << endl;
            }
        } 
    }
    input.close();
    output.close();

    return 0;
}
