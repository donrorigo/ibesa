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
    int j;
    int flag=atoi(argv[3]);
    double size;
    string file = argv[1], code, sequence;
    double lrcs;
    
    ifstream proteins;
    string amino_sequence;
    

    proteins.open(argv[2]);
    while(!proteins.eof())
    {
        getline(proteins, code, ' ');
        getline(proteins, amino_sequence, ' ');
        if(file == code) break;
        getline(proteins, amino_sequence);
    }
    proteins.close();

    std::vector<std::string> result;
    ofstream output("NORMALIZADO: " + file.substr(0,6));
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
                output << (stod(result[i])/0.35) << " ";
                break;
            case 2:
                if(flag == 1){
                    lrcs = stod(result[i]) / (amino_sequence.length()*3);
                    output << lrcs << endl;
                }
                else output << result[i] << endl;
            }
        } 
    }
    input.close();
    output.close();

    return 0;
}
