#include "ibesa.h"

/*
    Comentario: Tercer objetivo a cumplir en la resolución de la recombinación homóloga. 
    Autor: Rodrigo Puerto Pedrera
    Fecha: ***
    Fuente: Art. del tutor Miguel Ángel Vega Rodriguez - Multi-Objective Artificial Bee Colony for designing multiple genes
    encoding the same protein
*/

using namespace std;

map<char, vector<string> > amino_codons  = {
    {'A', {"GCA", "GCC", "GCG", "GCU"}},
    {'C', {"UGC", "UGU"}},
    {'D', {"GAC", "GAU"}},
    {'E', {"GAA", "GAG"}},
    {'F', {"UUC", "UUU"}},
    {'G', {"GGA", "GGC", "GGG", "GGU"}},
    {'H', {"CAC", "CAU"}},
    {'I', {"AUA", "AUC", "AUU"}},
    {'K', {"AAA", "AAG"}},
    {'L', {"CUA", "CUC", "CUG", "CUU", "UUA", "UUG"}},
    {'M', {"AUG"}},
    {'N', {"AAC", "AAU"}},
    {'P', {"CCA", "CCC", "CCG", "CCU"}},
    {'Q', {"CAA", "CAG"}},
    {'R', {"AGA", "AGG", "CGA", "CGC", "CGG", "CGU"}},
    {'S', {"AGC", "AGU", "UCA", "UCC", "UCG", "UCU"}},
    {'T', {"ACA", "ACC", "ACG", "ACU"}},
    {'V', {"GUA", "GUC", "GUG", "GUU"}},
    {'W', {"UGG"}},
    {'Y', {"UAC", "UAU"}}
    };

map <string, char> which_amino = {
    {"GCA", 'A'}, {"GCC", 'A'}, {"GCG", 'A'}, {"GCU", 'A'},
    {"UGC", 'C'}, {"UGU", 'C'},
    {"GAC", 'D'}, {"GAU", 'D'},
    {"GAA", 'E'}, {"GAG", 'E'},
    {"UUC", 'F'}, {"UUU", 'F'},
    {"GGA", 'G'}, {"GGC", 'G'}, {"GGG", 'G'}, {"GGU", 'G'},
    {"CAC", 'H'}, {"CAU", 'H'},
    {"AUA", 'I'}, {"AUC", 'I'}, {"AUU", 'I'},
    {"AAA", 'K'}, {"AAG", 'K'},
    {"CUA", 'L'}, {"CUC", 'L'}, {"CUG", 'L'}, {"CUU", 'L'}, {"UUA", 'L'}, {"UUG", 'L'},
    {"AUG", 'M'},
    {"AAC", 'N'}, {"AAU", 'N'},
    {"CCA", 'P'}, {"CCC", 'P'}, {"CCG", 'P'}, {"CCU", 'P'},
    {"CAA", 'Q'}, {"CAG", 'Q'},
    {"AGA", 'R'}, {"AGG", 'R'}, {"CGA", 'R'}, {"CGC", 'R'}, {"CGG", 'R'}, {"CGU", 'R'}, {"CAG", 'R'},
    {"AGC", 'S'}, {"AGU", 'S'}, {"UCA", 'S'}, {"UCC", 'S'}, {"UCG", 'S'}, {"UCU", 'S'},
    {"ACA", 'T'}, {"ACC", 'T'}, {"ACG", 'T'}, {"ACU", 'T'},
    {"GUA", 'V'}, {"GUC", 'V'}, {"GUG", 'V'}, {"GUU", 'V'},
    {"UGG", 'W'},
    {"UAC", 'Y'}, {"UAU", 'Y'}};

map<string, double> amino_weights = {
    /*aminoacid: A*/
    {"GCA", 0.390474084}, {"GCC", 0.532551795}, {"GCG", 0.136695421}, {"GCU", 1},
    /*aminoacid: C*/
    {"UGC", 0.404325033}, {"UGU", 1},
    /*aminoacid: D*/
    {"GAC", 0.703793889}, {"GAU", 1},
    /*aminoacid: E*/
    {"GAA", 1}, {"GAG", 0.315994266},
    /*aminoacid: F*/
    {"UUC", 1}, {"UUU", 0.942067628},
    /*aminoacid: G*/
    {"GGA", 0.177201478}, {"GGC", 0.229387027}, {"GGG", 0.118006882}, {"GGU", 1},
    /*aminoacid: H*/
    {"CAC", 0.761111111}, {"CAU", 1},
    /*aminoacid: I*/
    {"AUA", 0.26277856}, {"AUC", 0.683539061}, {"AUU", 1},
    /*aminoacid: K*/
    {"AAA", 0.846792801}, {"AAG", 1},
    /*aminoacid: L*/
    {"CUA", 0.310150799}, {"CUC", 0.093180284}, {"CUG", 0.21396954}, {"CUU", 0.240603196}, {"UUA", 0.64138345}, {"UUG", 1},
    /*aminoacid: M*/
    {"AUG", 1},
    /*aminoacid: N*/
    {"AAC", 1}, {"AAU", 0.872202532},
    /*aminoaci: P*/
    {"CCA", 1}, {"CCC", 0.184718349}, {"CCG", 0.11868377}, {"CCU", 0.510317903},
    /*aminoacid: Q*/
    {"CAA", 1}, {"CAG", 0.301447165},
    /*aminoacid: R*/
    {"AGA", 1}, {"AGG", 0.222301717}, {"CGA", 0.049979558}, {"CGC", 0.067252657}, {"CGG", 0.034955029}, {"CGU", 0.338000818},
    /*aminoacid: S*/
    {"AGC", 0.261645885}, {"AGU", 0.386334165}, {"UCA", 0.457157107}, {"UCC", 0.638703242}, {"UCG", 0.210673317}, {"UCU", 1},
    /*aminoacid: T*/
    {"ACA", 0.513350999}, {"ACC", 0.678760701}, {"ACG", 0.197513249}, {"ACU", 1},
    /*aminoacid: V*/ 
    {"GUA", 0.283953854}, {"GUC", 0.604002797}, {"GUG", 0.32337004}, {"GUU", 1},
    /*aminoacid: W*/
    {"UGG", 1},
    /*aminoacid: Y*/
    {"UAC", 1}, {"UAU", 0.810795614}};


void sorting_population()
/* ordena mediante fitness el vector pasado por referencia */
{
    std::sort(population.begin(), population.end(), [&](const single &x, const single &y) -> bool { return x.fitness > y.fitness; });
    return;
}

void create_superCAI(int numC, string sequence)
/* introduce en la población un elefante con CAI máximo */
{
    single result;
    int random_codon, i;
    string codon, id = "";

    for(int j=0; j < numC; ++j) /* numero de codones a generar */
        { 
            string CDS="";
            for(char aminoacid : sequence) /* bucle para generar el CDS */
            {               
                for(int i=0; i<amino_codons[aminoacid].size(); ++i)
                if(amino_weights[amino_codons[aminoacid][i]] == 1){
                    CDS += amino_codons[aminoacid][i];
                    break;
                }
            } 
            result.cds.push_back(CDS);
        }
    
    result.objetives.push_back(mCAI(result.cds).value);
    id += to_string(result.objetives[0]);
    result.objetives.push_back(mHD(result.cds).value);
    id += to_string(result.objetives[1]);
    result.objetives.push_back(mlrcs(result.cds).value);
    id += to_string(result.objetives[2]);
    result.id = id;
    result.fitness = 0;
    result.age = 0;
    result.gender = true; 
    result.lastmutation = -1;
    result.number = 0;
    population.push_back(result);

    return;
}

void init(int total, string amino_sequence, int CDSs, int machos)
/*  Inicialización de cada individuo de la población
    total = numero total de individuos
    file = archivo de partida 
    machos = numero de machos
    return -> lista population inicializada */
{
    int random_codon;
    int threads = omp_get_max_threads();
    string id;
     
           
    for(int i=1; i < 2*total; ++i) // numero de individuos de la población 
    {          
        single individuo;
        id = "";
        for(int j=0; j < CDSs; ++j) // numero de codones a generar 
        { 
            string CDS="";
            for(char aminoacid : amino_sequence) // bucle para generar el CDS 
            {                
                random_codon = rand_r(&random_vector[0]) % amino_codons[aminoacid].size(); // indice generado del codón a introducir 
                CDS += amino_codons[aminoacid][random_codon];
            } 
            individuo.cds.push_back(CDS);
        }
        
        individuo.objetives.push_back(mCAI(individuo.cds).value);
        id += to_string(individuo.objetives[0]);
        individuo.objetives.push_back(mHD(individuo.cds).value);
        id += to_string(individuo.objetives[1]);
        individuo.objetives.push_back(mlrcs(individuo.cds).value);
        id += to_string(individuo.objetives[2]);
        individuo.id = id;
        individuo.fitness = 0;
        individuo.age = 0;
        individuo.lastmutation = -1;
        individuo.gender = (i < machos) ? true : false; 
        individuo.number = i;
        population.push_back(individuo);
    }

    // for (int i = 0; i < total; ++i) population.push_back(population[1]);
       
    return;
}

void show_cdss(vector<string> CDSs)
/*  muestra por consola todos los CDS pasados por parámetro */
{
    for(int i=0; i<(int)CDSs.size(); ++i)
        std::cout << "CDS" << i << ": " << CDSs[i] << "\n\n";
    return;
}

void write_results()
/* escribe los elefantes con mejor fitness encontrados */
{

    ofstream fs(code + (".txt"));
    fs << fixed;
    fs << setprecision(15);
    int j, i=0;

    for(single solution : paretofront)
    {   
        //if(i==100) break;

        j=0;
        fs << "individual" << i << " mCAI=" << solution.objetives[0] << ",mHD=" <<  solution.objetives[1] 
        << ",lLCS=" <<  solution.objetives[2] << endl;
        for(string cds : solution.cds){
            fs << "CDS-" << j << "   " << cds << endl;
            j++;
        } 
        fs << "//" << endl;
        i++;
    }
    fs.close();
    return;
}

void show_population(int x)
/* muestra por consola toda la información de la población entera */
{
    fstream fs;
    
    fs.open(code + "_population.txt", std::fstream::out | std::fstream::app);
    int i=0, j;

    fs << "\nIteracción: " << x << endl;
    for(int n = 0; n < population.size(); ++n)
    {   
        if(n==population.size()/2) fs << endl;
        j=0;
        fs << "individual" << n << " mCAI=" << population[n].objetives[0] << ",mHD=" << population[n].objetives[1] << ",lLCS=" << population[n].objetives[2];
        fs << ",GENDER=" << ((population[n].gender) ? "male" : "female") << ",AGE:" << population[n].age <<",fitness:" << population[n].fitness << ", last mutation: " << population[n].lastmutation; 
        fs << ", Identificador de poblacion: " << population[n].number << endl;
        i++;
    }

    fs << "//" << endl;
    fs.close();
    return;
}

void export2utility()
/* exportación de los datos de debug al fichero utility.txt */
{
    ofstream fs(code + "_utility.txt");
    int machos, hembras;
    float m, h, o, g, opt;
    m=h=o=g=opt=0;
    
    fs << "Proteina: " << code << endl << "-------------" << endl;

    for(int i=0; i<howmany.size(); ++i)
    {
        machos = howmany[i].first;
        hembras = howmany[i].second;
        fs << "[!] Epoca "<< i << endl << endl;
        fs << "       Número de machos y hembras:" << endl;
        fs << "           - machos = " << machos << endl;
        fs << "           - hembras = " << hembras << endl;
        fs << "       Número de elefantes fallecidos por edad: " << oldest[i] << endl;
        fs << "       Número de «greedy mutations» que no sirven: " << nonutility[i] << endl;
        fs << "       Número de «optimum mutations» que no han servido: " << optimum_utility[i] << endl << endl;
        m+=howmany[i].first;
        h+=howmany[i].second;
        o+=oldest[i];
        g+=nonutility[i];
        opt+=optimum_utility[i];
    } 

    fs << "\n\n [?] PORCENTAJES DEL LANZAMIENTO: " << endl;
    fs << " [*] PORCENTAJE DE MACHOS: " << ((m/10000)*100) << endl;
    fs << " [*] PORCENTAJE DE HEMBRAS: " << ((h/10000)*100) << endl;
    fs << " [*] PORCENTAJE DE ELEFANTES FALLECIDOS: " << ((o/10000)*100) << endl;
    fs << " [*] PORCENTAJE DE 'GREEDY MUTATIONS' SIN UTILIDAD: " << ((g/h)*100) << endl;
    fs << " [*] PORCENTAJE DE 'OPTIMUM MUTATIONS' SIN UTILIDAD: " << ((opt/o)*100) << endl;
    fs << " [*] HIPERVOLUMEN OBTENIDO: " << endl;
    
    fs.close();
    return;
}

int three_mutations(int j, int id_th)
{
    double dcai, dmhd, dlrcs, dif;
    int iter;
    optimum_mutations[0](population[j], optimum_mutated[id_th][0], rand_r(&random_vector[id_th])% 5 + (GREEDYMUTATION-40), id_th, random_vector, auxiliar_cdss);
    optimum_mutations[1](population[j], optimum_mutated[id_th][1], rand_r(&random_vector[id_th])% 5 + (GREEDYMUTATION-40), id_th, random_vector, auxiliar_cdss);
    optimum_mutations[2](population[j], optimum_mutated[id_th][2], rand_r(&random_vector[id_th])% 5 + (GREEDYMUTATION-40), id_th, random_vector, auxiliar_cdss);
    dcai = std::abs(optimum_mutated[id_th][0].objetives[0] - population[j].objetives[0]);
    dmhd = std::abs(optimum_mutated[id_th][1].objetives[1] - population[j].objetives[1]);
    dlrcs = std::abs(optimum_mutated[id_th][2].objetives[2] - population[j].objetives[2]);
    dif = std::max({dcai, dmhd, dlrcs}, [](const std::double_t& s1, const std::double_t& s2) { return s1 < s2;});
    
    if(dif==dcai) return 0;
    else if(dif==dmhd) return 1;
    else return 2;
}

int main(int argc, char const *argv[])
/*  Argumentos: 
        1. número total de individuos de la población
        2. numero de épocas
        3. numero de individuos machos 
        4. codigo de la secuencia
        5. secuencia de aminoácidos
        6. numero de CDSs a generar
*/
{
    cout << "[*] Proteína " << code << " iniciada." << endl;

    /* definición de la variable de generaciones*/
    int i=0;

    /* guarda de los parámetros */
    int poblacion = atoi(argv[1]), epochs = atoi(argv[2]), total_cds = stoi(argv[6]);
    string aminoacids = argv[5];
    code = argv[4];

    /* definición de variables para el paralelismo */
    int id_th, hilos = omp_get_max_threads();

    /* inicialización del vector de llamadas a mutaciones */
    greedy_mutations = {lrcs_mutation, mhd_mutation, cai_mutation};
    optimum_mutations = {undue_cai_mutation, undue_mhd_mutation, undue_lrcs_mutation};

    /* reserva de memoria para el calculo del fitness */
    indicators.reserve(poblacion*2);
    bounds.reserve(3);
    for(int i=0; i<poblacion*2; ++i) indicators[i].reserve(poblacion*2);       
    
    /* inicialización de las estructuras de estadistica */
    volatile int machos, old, nonutil, optimum_util;
    for(int n=0; n<poblacion; ++n)
    {
        howmany.insert(std::make_pair(n, std::pair<int,int>(0,0)));
        oldest.insert(std::make_pair(n,0));
        nonutility.insert(std::make_pair(n,0));
        optimum_utility.insert(std::make_pair(n,0));
    }

    /* inicialización del vector de randoms */
    unsigned int seed = time(NULL); 
    for(int th=0; th<hilos; ++th) random_vector.push_back(seed+th);

    /* inicialización de la población */
    create_superCAI(total_cds, aminoacids);
    init(poblacion, aminoacids, total_cds, machos);
    for(int th=0; th<hilos; ++th) optimum_mutated.push_back({population[0],population[0],population[0]});        
    
    /* inicialización del vector auxiliar */
    for(int th = 0; th < hilos; ++th) auxiliar_cdss.push_back(population[0].cds);

    /* calculo del fitness y ordenacion de la primera poblacion */
    #pragma omp parallel
    compute_fitness(population, indicators, bounds);

    /* selección por elitismo (mayor fitness) */
    sorting_population();   
 
    /* inicio del algoritmo */
    #pragma omp parallel private(id_th) shared(machos, old, nonutil, optimum_util)
    {   
        id_th = omp_get_thread_num();
        while(i < epochs)   
        {
            #pragma omp single
            {
                show_population(i);
                cout << "[!] Generación: " << i  << endl;
                machos = old = nonutil = optimum_util = 0;
            }

            /* equilibrado de generos */
            #pragma omp for schedule(guided)
            for(int j=0; j<poblacion; ++j)  
            population[j].gender = (j<poblacion/2) ? false : true;

            /* mutaciones */
            #pragma omp for schedule(guided) reduction(+: nonutil)   
            for(int j=0; j<poblacion;j++)
            {
                if (!population[j].gender) {
                    greedy_mutations[rand_r(&random_vector[id_th])%3](population[j], population[j+poblacion], GREEDYMUTATION, id_th, random_vector, auxiliar_cdss);
                    
                    /* control de convergencia a una única solución */
                    if (population[j] == population[j+poblacion]) 
                    {
                        random_mutation(population[j], population[j+poblacion], 5*RANDOMMUTATION, id_th, random_vector, auxiliar_cdss); 
                        nonutil++;
                    }
                }else{
                    random_mutation(population[j], population[j+poblacion], RANDOMMUTATION, id_th, random_vector, auxiliar_cdss);
                    machos++;
                } 
                
                /* aumento de edad: si el nuevo elefante no domina al anterior se aumenta la edad del anterior */
                if(!(dominates(population[j+poblacion], population[j]) == 1)) population[j].age++;
            } 
            
            /* reset de los elefantes viejos */
            #pragma omp for schedule(guided) reduction(+: old, optimum_util)
            for(int j=0; j<2*poblacion; ++j)  /* mutaciones óptimas */
            { 
                if(population[j].age >= OLD) 
                {
                    old++;  
                    if(dominates(optimum_mutated[id_th][three_mutations(j, id_th)], population[j]) == 1) population[j]=optimum_mutated[id_th][0]; 
                    else
                    {
                        random_mutation(population[j], population[j], 10*RANDOMMUTATION, id_th, random_vector, auxiliar_cdss); 
                        optimum_util++;
                    }

                    /* reset de la edad */
                    population[j].age = 0;
                }
            }
            
            /* cálculo del fitness */
            compute_fitness(population, indicators, bounds);
            
            #pragma omp single
            {
                try
                {
                    sorting_population();   /* selección por elitismo (mayor fitness) */
                }
                catch(const std::exception& e)
                {
                    std::cerr << "[DANGER] Excepción en la ordenación población: "  << e.what() << '\n';
                    i=-1;
                }
            }

            #pragma omp single
            {
                /* guardado de la población mutada */
                std::copy(population.begin()+poblacion, population.end(), std::back_inserter(solutions)); 

                /* actualizado de vectores de utilidad */
                howmany[i].first = machos;
                howmany[i].second = poblacion - machos;
                oldest[i] = old;
                nonutility[i] = nonutil;
                optimum_utility[i] = optimum_util;

                /* aumento de generación */
                i++; 
            }
        }    
    }

    /* guardado de la última población */
    std::copy(population.begin(), population.end()-poblacion, std::back_inserter(solutions)); 

    /* creación del frente de pareto */
    for(i=0; i<solutions.size(); ++i)
    {
        dominated = false;

        /* control de elefante dominado */
        #pragma omp parallel for 
        for(int j=0; j<solutions.size(); ++j)
        {
            if(dominated) continue;
            if(dominates(solutions[j], solutions[i]) == 1) dominated = true; 
        } 

        /* control de no repetido */
        if(!dominated)
        {
            if(nonrepeat.find(solutions[i].id)==nonrepeat.end()) 
            {
                nonrepeat.insert(solutions[i].id);
                paretofront.push_back(solutions[i]);
            }
        }
    }

    /* escritura en fichero */
    write_results();
    export2utility();

    cout << "[*] Proteína " << code << " terminada.\n" << endl; 

    return 0;
}