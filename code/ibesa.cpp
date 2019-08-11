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


void sorting_population(vector<single> &vector_pop)
/* ordena mediante fitness el vector pasado por referencia */
{
    std::sort(vector_pop.begin(), vector_pop.end(), [&](const single &x, const single &y) -> bool {
        return x.fitness > y.fitness;
    });

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
    string id;
    create_superCAI(CDSs,amino_sequence);

    for(int i=1; i < total*2; ++i) // numero de individuos de la población 
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
        individuo.gender = (i < machos) ? true : false; 
        population.push_back(individuo);
        
    }

    
    return;
}

void stabilize_population()
/* permite la no estancación del algoritmo cambiando el sexo de los elefantes */
{

    for(size_t e=0; e<population.size(); ++e){
        
        if(population[e].objetives[0] > 0.65 && population[e].objetives[0] < 0.9 && fitness_vector.find((population[e].id))==fitness_vector.end()){
            //solutions.push_back(population[e]);
            //fitness_vector.insert(value);
            cout << "entré! CAI" << "objetivos: " << population[e].id << endl;
        } 
        else if(population[e].objetives[1] > 0.22 && fitness_vector.find((population[e].id))==fitness_vector.end())
        {
            //fitness_vector.insert(value);
            //solutions.push_back(population[e]);
            cout << "entré! MHD" << "objetivos: " << population[e].id << endl;
        }
        else if(population[e].objetives[0] < 0.10 && fitness_vector.find((population[e].id))==fitness_vector.end())
        {
            //fitness_vector.insert(value);
            //solutions.push_back(population[e]);
            cout << "entré! LRCS" << "objetivos: " << population[e].id << endl;
        } 
    }
}

void show_cdss(vector<string> CDSs)
/*  Muestra por consola todos los CDS pasados por parámetro */
{
    for(int i=0; i<(int)CDSs.size(); ++i)
        std::cout << "CDS" << i << ": " << CDSs[i] << "\n\n";
    return;
}

void write_results(string code)
/* escribe los elefantes con mejor fitness encontrados */
{

    ofstream fs(code + (".txt"));
    fs << fixed;
    fs << setprecision(15);
    int j, i=0;

    for(single solution : solutions)
    {   
        if(i==100) break;

        j=0;
        fs << "individual" << i << " mCAI=" << solution.objetives[0] << ",mHD=" <<  solution.objetives[1] 
        << ",lLCS=" <<  solution.objetives[2] << ",Fitness=" << solution.fitness << ",Gender=" << (solution.gender ? "Male" : "Female") << endl;
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

void show_single(single a)
/* muestra por consola toda la información del elefante pasado por parámetro */
{
    int j= 0;
    std::cout << " mCAI=" << a.objetives[0] << ",mHD=" << a.objetives[1] << ",lLCS=" << a.objetives[2];
    std::cout << ",GENDER=" << ((a.gender) ? "male" : "female") << a.age <<",fitness:" << a.fitness << endl;
    for(string cds : a.cds){
        std::cout << "CDS-" << j << "   " << cds << endl;
        j++;
    } 
}

void show_population()
/* muestra por consola toda la información de la población entera */
{

    int i=0, j;
    single result;

    for(int n = 0; n<population.size(); ++n)
    {   
        result = population[n];
        j=0;
        std::cout << "individual" << i << " mCAI=" << result.objetives[0] << ",mHD=" << result.objetives[1] << ",lLCS=" << result.objetives[2];
        std::cout << ",GENDER=" << ((result.gender) ? "male" : "female") << ",AGE:" << result.age <<",fitness:" << result.fitness << endl;
        for(string cds : result.cds){
            std::cout << "CDS-" << j << "   " << cds << endl;
            j++;
        } 
        std::cout << "//" << endl;
        i++;
    }

    return;
}

int main(int argc, char const *argv[])
/*  Argumentos: 
        1. número total de individuos de la población
        2. numero de épocas
        3. numero de individuos machos 
        4. codigo de la secuencia
        5. secuencia de aminoácidos
        6. numero de CDSs a generar*/
{
    /* definición de variables e inicialización */
    int j, i=0;
    int poblacion = atoi(argv[1]), epochs = atoi(argv[2]), machos = atoi(argv[3]);
    greedy_mutations = {cai_mutation, mhd_mutation, lrcs_mutation, undue_cai_mutation};
    unsigned int seed = time(NULL); 
    int id_th, hilos = omp_get_max_threads();
    int total_cds = stoi(argv[6]);
    string aminoacids = argv[5];

    /* reserva de memoria */
    population.reserve(poblacion*2);
    bounds.reserve(3);
    indicators.reserve(poblacion*2);
    for(int i=0; i<poblacion*2; ++i) indicators[i].reserve(poblacion*2);       
    random_vector.reserve(hilos);
    int len = 3 * aminoacids.length();
    for(int th=0; th<hilos; ++th) random_vector[th] = seed+th;

    /* inicialización de la población */
    init(poblacion, aminoacids, total_cds, machos);

    auxiliar_cdss.resize(hilos);
    for(int th = 0; th < hilos; ++th) auxiliar_cdss[th] = population[0].cds;

    #pragma omp parallel private(id_th)
    {   
        id_th = omp_get_thread_num();     

        while(i < epochs)
        {
            #pragma omp single
            cout << "Epoca: " << i << endl;
            
            #pragma omp single
            cout << "[!] Mutaciones." << endl;
            
            #pragma omp for schedule(guided)  
            for(j=0; j<poblacion;j++)
            {
                // (population[j].gender) ? random_mutation(population[j], population[j+poblacion], 60, id_th, random_vector, auxiliar_cdss) : 
                // greedy_mutations[rand_r(&random_vector[id_th])%2](population[j], population[j+poblacion], 83, id_th, random_vector, auxiliar_cdss);
                cai_mutation(population[j], population[j+poblacion], 83, id_th, random_vector, auxiliar_cdss);
                if(dominates(population[j+poblacion], population[j]) == -1) population[j].age = 0;
                else population[j].age++;
            } 

            #pragma omp single
            cout << "[!] Control de elefantes mayores." << endl;
            
            #pragma omp single
            for(j=0; j<poblacion; ++j){
                if(population[j].age == OLD){
                    if(fitness_vector.find((population[j].id))==fitness_vector.end())
                    {
                        cout << "Entre! Elefante mayor!" << "objetivos: " << population[j].id << endl;
                        cout << "Elefante tiene estos objetivos: " << population[j].objetives.size() << endl;
                        //solutions.push_back(population[j]);
                        //fitness_vector.insert(value);
                    }
                    
                } 
            }

            #pragma omp single
            cout << "[!] Criba de elefantes mayores-" << endl;

            #pragma omp single
            cout << "[!] Tamaño del vector de objetivos: " << fitness_vector.size() << endl;
            
            #pragma omp for schedule(guided)  
            for(j=0; j<poblacion; ++j){
                if(population[j].age == OLD){
                    random_mutation(population[j], population[j], 90, id_th, random_vector, auxiliar_cdss);
                    population[j].gender = true;
                }
            }
            
            #pragma omp single
            cout << "[!] Fitness." << endl;
            
            compute_fitness(population, indicators, bounds);
            #pragma omp single
            {
                cout << "[!] Ordenación" << endl;
            
                sorting_population(population);
                stabilize_population();
                cout << fitness_vector.size() << endl;
            
                i++;     
            }
        }
    
        #pragma omp single
        {
            population.clear();
            indicators.resize(solutions.size());
        }

        #pragma omp for schedule(guided)
        for(int i=0; i<solutions.size(); ++i) indicators[i].resize(solutions.size());   
        compute_fitness(solutions, indicators, bounds);    
        
    }

    sorting_population(solutions);
    write_results(argv[4]);
    for(string sda : fitness_vector) cout << sda << endl;
    
    return 0;
}