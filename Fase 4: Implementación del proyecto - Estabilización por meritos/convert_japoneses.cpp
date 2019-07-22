#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstring>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char const *argv[])
{
    int j;
    bool flag;
    string file = argv[1], code, cai, mhd, lrcs;
    ofstream output("RESULTADOS: " + file.substr(0,6));
    ifstream input;
    
    input.open(file);
    while(!input.eof())
    {
        getline(input, code);
        if(code[0] != 'i') continue;
        j=0;
        flag = false;
        cai = mhd = lrcs = "";
        for(int i =0; i<code.length(); i++){
            switch (j)
            {
            case 0:
                if(code[i]==',') j++;
                else cai +=  code[i];
                break;
            
            case 1:
                if(code[i]==',') j++;
                else mhd += code[i];
                break;
            
            case 2:
                if(code[i]==',') j++;
                else lrcs += code[i];
                break;
            }
        } 


      
        for(int i =0; i<cai.length(); i++){
            if(cai[i] == '='){
                flag = true;
                i++;
            } 
            if(flag) output << cai[i]; 
            
        } 
        
        flag = false;
        output << " ";

        for(int i =0; i<mhd.length(); i++){
            if(mhd[i] == '='){
                flag = true;
                i++;
            } 
            if(flag) output << mhd[i]; 
            
        } 

        flag = false;
        output << " ";
         
        for(int i=0; i<lrcs.length()-2; i++){
            if(isdigit(lrcs[i])) output <<   lrcs[i];
        }
        
        output << " \n";
    }

    input.close();
    output.close();

    return 0;
}
