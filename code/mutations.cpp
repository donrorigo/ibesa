#include "mutations.h"

/*
    Comentario: Tercer objetivo a cumplir en la resolución de la recombinación homóloga. 
    Autor: Rodrigo Puerto Pedrera
    Fecha: mar jun 18, 2019
    Fuente: Art. del tutor Miguel Ángel Vega Rodriguez - Multi-objective protein encoding: Redefinition of the problem,
    new problem-aware operators, and approach based on Variable Neighborhood Search.
*/

void random_mutation(single &a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss)
/*  compute a new single */
{
    string codon, random_codon, new_cds = "", id = "";
    string CDS;
    int size, c=0;
    
    //printf("[*] Hilo %d: Mutación random.\n", th);

    for(int c=0; c<a.cds.size(); ++c)
    {
        CDS = a.cds[c]; 
        for(int i = 0; i<(int)CDS.size(); i+=3)
        {
            codon = CDS.substr(i,3);
            size = amino_codons[which_amino[codon]].size();

            if(rand_r(&random_vector[th]) % 100 < Pm && size > 1)
            {
                do{
                    random_codon = amino_codons[which_amino[codon]][rand_r(&random_vector[th])%size];
                    if(random_codon!=codon) new_cds += random_codon;
                }while(random_codon==codon);
            }else new_cds += codon;
        }

        result.cds[c] = (new_cds);
        new_cds = "";
    }

    result.objetives[0] = (mCAI(result.cds).value);
    id += to_string(result.objetives[0]);
    result.objetives[1] = (mHD(result.cds).value);
    id += to_string(result.objetives[1]);
    result.objetives[2] = (mlrcs(result.cds).value);
    id += to_string(result.objetives[2]);
    result.id = id;
    result.fitness = 0;
    result.age = a.age;
    result.gender = a.gender;
    result.number = a.number;
    result.lastmutation = 0;
    
    return;
}

void cai_mutation(single & a, single & result, double Pm, int th, vector<unsigned int> & random_vector, vector<vector<string> > & auxiliar_cdss)
/*  primera mutación avariciosa */
{
    string codon, random_codon, new_cds = "", CDS = mCAI(a.cds).cds1, id = "";

    //printf("[*] Hilo %d: Mutación avariciosa CAI.\n", th);
    try
    {
        for(int i = 0; i<(int)CDS.size(); i+=3)
        {
            codon = CDS.substr(i,3);   
            if(rand_r(&random_vector[th]) % 100 < Pm && amino_weights[codon] != 1)  
            {
                do{
                    random_codon = amino_codons[which_amino[codon]][rand_r(&random_vector[th])%amino_codons[which_amino[codon]].size()];
                    if(amino_weights[random_codon] > amino_weights[codon])  new_cds += random_codon;                
                }while(amino_weights[random_codon] < amino_weights[codon] || random_codon == codon);
            }else{
                new_cds += codon;
            }      
        }

        update_CDSs(a.cds, new_cds, CDS, auxiliar_cdss[th]);
        update_vector(auxiliar_cdss[th], result.cds);
        result.objetives[0] = (mCAI(result.cds).value);
        id += to_string(result.objetives[0]);
        result.objetives[1] = (mHD(result.cds).value);
        id += to_string(result.objetives[1]);
        result.objetives[2] = (mlrcs(result.cds).value);
        id += to_string(result.objetives[2]);
        result.id = id;
        result.fitness = 0;
        result.age = a.age;
        result.gender = a.gender;
        result.number = a.number;
        result.lastmutation = 1;

        
    }
    catch(const std::exception& e)
    {
        std::cerr << "[DANGER] Mutación avariciosa CAI ha provocado la siguiente excepción: " << e.what() << '\n';
    }
    


    return;
}

void mhd_mutation(single & a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss)
{
    result.cds = a.cds;
    aim aim_v = mHD(a.cds);   
    string codon, best_codon, new_CDS1, id = "";
    aim curr_HD, curr_mHD=aim_v, best_HD, best_mHD, new_HD, new_mHD;
    best_HD.value=-1;
    best_mHD.value=-1;
    try
    {
        for(int i = 0; i<(int)aim_v.cds1.size(); i+=3) 
        {
            codon = aim_v.cds1.substr(i,3);
            if(rand_r(&random_vector[th]) % 100 < Pm &&  amino_codons[which_amino[codon]].size() > 1)
            {
                curr_HD.value = HD(aim_v.cds1,aim_v.cds2);

                for(string new_codon : amino_codons[which_amino[codon]])
                {
                    if(new_codon != codon)
                    {
                        new_CDS1 = change_CDS(aim_v.cds1, new_codon, i);
                        new_HD.value = HD(new_CDS1, aim_v.cds2);
                        update_CDSs(a.cds, new_CDS1, aim_v.cds1, auxiliar_cdss[th]);
                        new_mHD = mHD(auxiliar_cdss[th]);
                        if(new_mHD.value > curr_mHD.value && new_mHD.value > best_mHD.value)
                        {
                            best_mHD = new_mHD;
                            best_codon = new_codon;
                        }else if (best_mHD.value == -1 && new_mHD.value == curr_mHD.value && new_HD.value > curr_HD.value && new_HD.value > best_HD.value)
                        {
                            best_HD = new_HD;
                            best_codon = new_codon;
                        }
                    }
                }

                if(best_mHD.value != -1 || best_HD.value != -1){
                    new_CDS1 = change_CDS(aim_v.cds1, best_codon, i);
                    update_CDSs(a.cds, new_CDS1, aim_v.cds1, auxiliar_cdss[th]);
                    update_vector(auxiliar_cdss[th], result.cds); 
                }

                best_mHD.value = -1;
                best_HD.value = -1;
            }

        }

        result.objetives[0] = (mCAI(result.cds).value);
        id += to_string(result.objetives[0]);
        result.objetives[1] = (mHD(result.cds).value);               
        id += to_string(result.objetives[1]);
        result.objetives[2] = (mlrcs(result.cds).value);
        id += to_string(result.objetives[2]);
        result.id = id;
    }
    catch(const std::exception& e)
    {
        std::cerr << "[DANGER] Mutación avariciosa MHD ha provocado la siguiente excepción: " << e.what() << '\n';
    }
    
    //printf("[*] Hilo %d: Mutación avariciosa MHD.\n", th);

    result.fitness = 0;
    result.age = a.age;
    result.gender = a.gender;
    result.number = a.number;
    result.lastmutation = 2;

    return;
}

void lrcs_mutation(single & a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss)
{
    aim curr_lrcs = mlrcs(a.cds), new_lrcs;
    string codon, random_codon, new_cds, cds1, cds2, cds3, CDS, id = "";
    result.cds = a.cds;
    int ret, index;

    try
    {
        for(char &character: curr_lrcs.cds1) cds1 += toupper(character);
        for(char &character: curr_lrcs.cds2) cds2 += toupper(character);  
        for(char &character: curr_lrcs.cds3) cds3 += toupper(character);   

        for(int i = 0; i<(int)curr_lrcs.cds1.size(); i+=3)
        {   
            /* cálculo del índice donde se encuentra, el CDS y el codon correspondientes */
            ret = ((curr_lrcs.index!=0) ? a.cds[0].length()%curr_lrcs.index : 0)%3;
            CDS = (curr_lrcs.index>=a.cds[0].length()) ? cds3 : cds2;
            index = (curr_lrcs.index>=a.cds[0].length()) ? curr_lrcs.index-(a.cds[0].length()) : curr_lrcs.index;
            codon = (ret!=0) ? CDS.substr((index+i-ret),3) : CDS.substr((index+i),3);
    
            if(rand_r(&random_vector[th]) % 100 < Pm && amino_codons[which_amino[codon]].size() > 1)
            {
                for(string random_codon : amino_codons[which_amino[codon]])
                {
                    if(random_codon != codon)
                    {
                        new_cds = ""; 
                        new_cds = change_CDS(CDS, random_codon, ((ret!=0) ? (index+i-ret): (index+i)));
                        update_CDSs(a.cds, new_cds, CDS, auxiliar_cdss[th]);
                        new_lrcs = mlrcs(auxiliar_cdss[th]);
                        if(new_lrcs.value < curr_lrcs.value)
                        {
                            update_vector(auxiliar_cdss[th], result.cds);
                            curr_lrcs.value = new_lrcs.value;    
                        }
                    }
                }
            }

        }

        result.objetives[0] = (mCAI(result.cds).value);
        id += to_string(result.objetives[0]);
        result.objetives[1] = (mHD(result.cds).value);
        id += to_string(result.objetives[1]);
        result.objetives[2] = curr_lrcs.value;
        id += to_string(result.objetives[2]);
        result.id = id;
    }
    catch(const std::exception& e)
    {
        std::cerr << "[DANGER] Mutación avariciosa LRCS ha provocado la siguiente excepción: " << e.what() << '\n';
    }
    
    //printf("[*] Hilo %d: Mutación avariciosa LRCS.\n", th);
    
    result.fitness = 0;
    result.age = a.age;
    result.gender = a.gender;
    result.number = a.number;
    result.lastmutation = 3;
    
    return;
}

                    

void undue_cai_mutation(single &a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss)
{
    string codon, random_codon, new_cds = "", id = "", new_CDS1;
    string CDS;
    int size, c=0;
    result.cds = a.cds;
    aim new_mHD, new_lrcs;

    try
    {
        /* for each cds of the protein*/
        for(int c=0; c<a.cds.size(); ++c)
        {
            CDS = a.cds[c]; 

            /* for each codon in the CDS */
            for(int i = 0; i<(int)CDS.size(); i+=3)
            {
                codon = CDS.substr(i,3);
                size = amino_codons[which_amino[codon]].size();

                /* if the P() < Pm and the number of synonyms is greater than itself */
                if(rand_r(&random_vector[th]) % 100 < Pm && amino_weights[codon] != 1)
                {
                    /* get the codon with CAI=1*/
                    do{
                        random_codon = amino_codons[which_amino[codon]][rand_r(&random_vector[th])%size];
                        if(amino_weights[random_codon] == 1)  new_cds += random_codon;                
                    }while(amino_weights[random_codon] != 1);

                    /* generate new CDS with the new codon */
                    new_CDS1 = change_CDS(CDS, random_codon, i);
                    update_CDSs(result.cds, new_CDS1, CDS, auxiliar_cdss[th]);

                    /* compute the new metrics */
                    new_mHD = mHD(auxiliar_cdss[th]);
                    new_lrcs = mlrcs(auxiliar_cdss[th]);

                    /* if any metric is better than the old one, we update the old CDS with the new one */
                    if(new_mHD.value >= a.objetives[1] && new_lrcs.value <= a.objetives[2])
                    {
                        update_vector(auxiliar_cdss[th], result.cds);
                        a.objetives[1] = new_mHD.value;
                        a.objetives[2] =  new_lrcs.value;
                    } 
                }
            }
            
        }

        /* update all the attributes of the elephant */
        result.objetives[0] = (mCAI(result.cds).value);
        id += to_string(result.objetives[0]);
        result.objetives[1] = (mHD(result.cds).value);               
        id += to_string(result.objetives[1]);
        result.objetives[2] = (mlrcs(result.cds).value);
        id += to_string(result.objetives[2]);
        result.id = id;
        result.fitness = 0;
        result.age = a.age;
        result.gender = a.gender;
        result.number = a.number;
        result.lastmutation = 4;
        return;
    }
    catch(const std::exception& e)
    {
        std::cerr << "[DANGER] Mutación óptima CAI ha provocado la siguiente excepción: " << e.what() << '\n';
    }
}

void undue_mhd_mutation(single &a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss)
/*  segunda mutación avariciosa */
{
    result.cds = a.cds;
    aim aim_v = mHD(a.cds);   
    string codon, best_codon, new_CDS1, id = "";
    aim curr_HD, curr_mHD=aim_v, best_HD, best_mHD, new_HD, new_mHD;
    best_HD.value=-1;
    best_mHD.value=-1;

    try
    {
        /* for each codon of the CDS returned by mHD function */
        for(int i = 0; i<(int)aim_v.cds1.size(); i+=3) 
        {
            codon = aim_v.cds1.substr(i,3);

            /* if P(t) & the codon's CDS has severals codons */
            if(rand_r(&random_vector[th]) % 100 < Pm &&  amino_codons[which_amino[codon]].size() > 1) 
            {
                curr_HD.value = HD(aim_v.cds1, aim_v.cds2);

                /* for each new codon of the mHDcodon's CDS */
                for(string new_codon : amino_codons[which_amino[codon]]) 
                {
                    if(new_codon != codon) 
                    {
                        /* new calculate with the new CDS constructed with the new codon*/
                        new_CDS1 = change_CDS(aim_v.cds1, new_codon, i);
                        new_HD.value = HD(new_CDS1, aim_v.cds2);
                        update_CDSs(a.cds, new_CDS1, aim_v.cds1, auxiliar_cdss[th]);
                        new_mHD = mHD(auxiliar_cdss[th]); 

                        /* if this improve mHD value, we update the mHD value */
                        if(new_mHD.value > curr_mHD.value && new_mHD.value > best_mHD.value)
                        {
                            best_mHD = new_mHD;
                            best_codon = new_codon;

                        /* else if this improve the HD value, we update the HD value */
                        }else if (best_mHD.value == -1 && new_mHD.value == curr_mHD.value && new_HD.value > curr_HD.value && new_HD.value > best_HD.value)
                        {
                            best_HD = new_HD;
                            best_codon = new_codon;
                        }
                    }
                }

                /* if some objective improves, we update the CDS of the elephant and the currently values*/
                if(best_mHD.value != -1 || best_HD.value != -1){
                    new_CDS1 = change_CDS(aim_v.cds1, best_codon, i);
                    update_CDSs(a.cds, new_CDS1, aim_v.cds1, auxiliar_cdss[th]);
                    update_vector(auxiliar_cdss[th], result.cds); 
                    aim_v = (best_mHD.value != -1) ? best_mHD : best_HD;
                    i = 0;
                }

                /* we reset the values for the next iteraction */
                best_mHD.value = -1;
                best_HD.value = -1;
            }

        }

        /* when we finish, we recalculate all the objectives for the elephant */
        result.objetives[0] = (mCAI(result.cds).value);
        id += to_string(result.objetives[0]);
        result.objetives[1] = curr_mHD.value;
        id += to_string(result.objetives[1]);
        result.objetives[2] = (mlrcs(result.cds).value);
        id += to_string(result.objetives[2]);
        result.id = id;
        result.fitness = 0;
        result.age = a.age;  
        result.gender = a.gender;
        result.number = a.number;
        result.lastmutation = 5;

        return;
    }
    catch(const std::exception& e)
    {
        std::cerr << "[DANGER] Mutación óptima MHD ha provocado la siguiente excepción: " << e.what() << '\n';
    }
}

void undue_lrcs_mutation(single &a, single &result, double Pm, int th, vector<unsigned int> &random_vector, vector<vector<string> > &auxiliar_cdss)
/*  tercera mutación avariciosa */
{
    aim curr_lrcs = mlrcs(a.cds), new_lrcs;
    string codon, random_codon, new_cds, cds1, cds2, cds3, CDS, id = "";
    result.cds = a.cds;
    int ret, index;
    
    try
    {
        /* convert the strings into uppercase characters for the correct performance of lcs and lrs algorithms */
        for(char &character: curr_lrcs.cds1) cds1 += toupper(character);
        for(char &character: curr_lrcs.cds2) cds2 += toupper(character);  
        for(char &character: curr_lrcs.cds3) cds3 += toupper(character);   

        for(int i = 0; i<(int)curr_lrcs.cds1.size(); i+=3) /* for each codon */
        {   
            /* which CDS of the elephant we choose */
            ret = ((curr_lrcs.index!=0) ? a.cds[0].length()%curr_lrcs.index : 0)%3;
            CDS = (curr_lrcs.index>=a.cds[0].length()) ? cds3 : cds2;
            index = (curr_lrcs.index>=a.cds[0].length()) ? curr_lrcs.index-(a.cds[0].length()) : curr_lrcs.index;
            codon = (ret!=0) ? CDS.substr((index+i-ret),3) : CDS.substr((index+i),3);

            if(amino_codons[which_amino[codon]].size() > 1) /* if the codon has synonyms codons*/
            {
                for(string random_codon : amino_codons[which_amino[codon]]) /* for each synonym codon */
                {
                    if(random_codon != codon) /* if its not the same we studied */
                    {
                        /* the new constructed CDS with the new codon */
                        new_cds = change_CDS(CDS, random_codon, ((ret!=0) ? (index+i-ret): (index+i)));
                        update_CDSs(a.cds, new_cds, CDS, auxiliar_cdss[th]);

                        /* new lrcs value */
                        new_lrcs = mlrcs(auxiliar_cdss[th]);

                        /* if the new lrcs value its better than the currently, we update it */
                        if(new_lrcs.value < curr_lrcs.value)
                        {
                            update_vector(auxiliar_cdss[th], result.cds);

                            /* we update the studied lrcs' CDSs for improve the improved */
                            curr_lrcs = new_lrcs;

                            /* reset the cycle */
                            i=0;   
                        }
                    }
                }
            }

        }
        
        /* we reset the values for the objectives */
        result.objetives[0] = (mCAI(result.cds).value);
        id += to_string(result.objetives[0]);
        result.objetives[1] = (mHD(result.cds).value);
        id += to_string(result.objetives[1]);
        result.objetives[2] = curr_lrcs.value;
        id += to_string(result.objetives[2]);
        result.id = id;
        result.fitness = 0;
        result.age = a.age;
        result.gender = a.gender;
        result.number = a.number;
        result.lastmutation = 6;
        
        return;
    }
    catch(const std::exception& e)
    {
        std::cerr << "[DANGER] Mutación óptima LRCS ha provocado la siguiente excepción: " << e.what() << '\n';
    }
}

string change_CDS(string cds, string new_codon, int index)
/* update the cds passed by parameter, changing the currently codon by the new one */
{
    string result = "";
    for(int i = 0; i<(int)cds.size(); i+=3) result += (i==index) ? new_codon : cds.substr(i,3);
    return result;    
}

void update_CDSs(vector<string> CDSs, string new_cds, string curr_cds, vector<string> & vector_aux)
/* update the vector passed by parameter, changing the currently cds by the new one */
{
    for(int i=0; i<CDSs.size(); ++i) vector_aux[i] = ((CDSs[i] == curr_cds) ? new_cds : CDSs[i]);  
    
    return;
}

void update_vector(vector<string> &src, vector<string> &dst)
/* update in memory the two vectors for the correct permormance in parallel compute */
{
    for(int i=0; i<src.size(); ++i) 
        for(int j=0; j<src[i].size(); ++j)
            dst[i][j] = src[i][j];
        
    return;
}