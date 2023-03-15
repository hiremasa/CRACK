#include <getopt.h> //options
#include <stdlib.h> //set handler
#include <iostream>
#include <vector>
#include <sys/types.h>
#include <fstream>
#include <stdlib.h> // rand

#include "rack_1.h"
#include "rack_2.h"
#include "synthetic_data.h"


void outofmemory(){
    fprintf(stderr, "Oh noes! Out of memory\n");
    exit(1);
}

void printOutput(pns::Result r, std::string file){
    std::pair< std::string, std::vector< pns::dot_file* > > out;
    std::string comp = "Compression: " + std::to_string(r.finalScore) + " / " + std::to_string(r.initialScore) + " = " + std::to_string(r.getCompression());
    out.first = comp;
    out.second = r.files;
    std::vector< pns::dot_file* > files = out.second;
    std::ofstream myfile;
    myfile.open(file.c_str());
    myfile << "digraph G {\n";
    myfile << "\tlabel = \"" << out.first << "\";\n";
    for(int i = 0; i < files.size(); i++){
        for(int j = 0; j < files[i]->size(); j++){
            if(myfile.is_open())
                myfile << (*(files[i]))[j];
        }
    }
    myfile << "}\n";
    myfile.close();
    
}

void deleteDotFiles(pns::Result r){
    for(int i = 0; i < r.files.size(); i++){
        delete r.files[i];
    }
}

bool dirExists(const char* file){
    std::ofstream myfile;
    myfile.open(file);
    if(myfile.is_open()){
        myfile.close();
        std::remove(file);
        return true;
    }else{
        fprintf(stderr, "Oh noes! Directory %s could not be opened or found!\n", file);
        return false;
    }
}

std::vector<int> randomSample(int max, int rows){
    std::vector<int> sample;
    for(int i = 0; i < rows; i++){
        int newOne = rand() % max;
        sample.push_back(newOne);
    }
    return sample;
}

/**
 * returns: 0 for --, 1 for -> and -1 for <-
 **/
std::pair<int,double> runCracko(Data *d, int precision, const char * append, const char *outfolder, bool bagging, bool norm, int argc, char** argv){
    pns::Result result3 = ERGO::runGreedyAlgorithm(d, precision, 2);
    pns::Result result4 = ERGO::runGreedyAlgorithm(d, precision, 3);
    double Lx = result3.initialScore;
    double Ly = result4.initialScore;
    double Lxy = result3.finalScore;
    double Lyx = result4.finalScore;
    double deltaYX = (Ly + Lxy) / (Lx + Ly);
    double deltaXY = (Lx + Lyx) / (Lx + Ly);
    double pval = std::pow(2.0,(-(std::abs((Lx + Lyx) - (Ly + Lxy))/2.0)));
    fprintf(stdout, "Lx|y:\t%f\nLy|x:\t%f\n", Lxy, Lyx);
    if (outfolder != NULL && !append && !bagging) {
        printOutput(result3, std::string(outfolder) + "Lxy.dot");
        printOutput(result4, std::string(outfolder) + "Lyx.dot");
    }
    double runtime = result3.runtime + result4.runtime;
    double epsilon = std::abs(deltaXY - deltaYX);
    //epsilon = pval;
    std::ofstream myfile;
    std::string output = "deltaX->Y: " + std::to_string(deltaXY);
    output += "\ndeltaY->X: " + std::to_string(deltaYX);
    std::string direction;
    std::pair<int,double> out;
    out.second = epsilon;
    //if(d->rows == 2)
    //    undecided = norm ? (deltaXY > 0.99 && deltaYX > 0.99) : false;
    if(epsilon == 0.0){ //pval >= 0.05 ||
        output += "\nResult: X--Y with epsilon = " + std::to_string(epsilon);
        direction = "--";
        out.first = 0;
    }else if(deltaXY < deltaYX){
        output += "\nResult: X->Y with epsilon = " + std::to_string(epsilon);
        direction = "->";
        out.first = 1;
    }else{
        output += "\nResult: X<-Y with epsilon = " + std::to_string(epsilon);
        direction = "<-";
        out.first = -1;
    }
    output += "\nRuntime: " +std::to_string(runtime) + " second\n";
    std::cout << "---------------------------\n" << output;
    if(outfolder != NULL && !append && !bagging){
        myfile.open((std::string(outfolder) + "summary.txt").c_str());
        if(myfile.is_open()){
            myfile << "Execution:\t";
            for(int i = 0; i < argc; i++){
                myfile << std::string(argv[i]) + " ";
            }
            myfile << "\n" + output;
        }
        myfile.close();
    }
    if(append && !bagging){
        myfile.open((std::string(outfolder) + "results.tab").c_str(), std::ios_base::app);
        if(myfile.is_open()){
            myfile << std::string(append) + "\t" + std::to_string(epsilon) + "\t" + direction + "\t" + std::to_string(runtime) + "\n";
        }
        myfile.close();
    }
    // delete dot files
    deleteDotFiles(result3);
    deleteDotFiles(result4);
    return out;
}

std::pair<int,double> runCracker2(Data *d, int precision, const char * append, const char *outfolder, bool bagging, bool norm, int argc, char** argv){
    d->nci = true;
    pns::Result result = ERGO::runGreedyAlgorithm(d, precision, 1);
    deleteDotFiles(result);
    std::pair<int,double> out;
    return out;
}


/**
 * returns: 0 for --, 1 for -> and -1 for <-
 **/
std::pair<int,double> runCracker(Data *d, int precision, const char * append, const char *outfolder, bool bagging, bool norm, int argc, char** argv){
    d->nci = true;
    pns::Result result2 = ERGO::runGreedyAlgorithm(d, precision, 2);
    pns::Result result3 = ERGO::runGreedyAlgorithm(d, precision, 3);
    double Lx = result2.initialScore;
    double Ly = result3.initialScore;
    double Lxy = result2.finalScore;
    double Lyx = result3.finalScore;
    fprintf(stdout, "Lx:\t%f\nLx|y:\t%f\nLy:\t%f\nLy|x:\t%f\n", Lx, Lxy, Ly, Lyx);
    double deltaXY = Lyx;
    if(norm)
        deltaXY /= Ly;
    double deltaYX = Lxy;
    if(norm)
        deltaYX /= Lx;
    if (outfolder != NULL && !append && !bagging) {
        printOutput(result2, std::string(outfolder) + "Lxy.dot");
        printOutput(result3, std::string(outfolder) + "Lyx.dot");
    }
    double runtime = result2.runtime + result3.runtime;
    double epsilon = std::abs(deltaXY - deltaYX);
    std::ofstream myfile;
    std::string output = "deltaX->Y: " + std::to_string(deltaXY);
    output += "\ndeltaY->X: " + std::to_string(deltaYX);
    std::string direction;
    std::pair<int,double> out;
    out.second = epsilon;
    bool undecided = false;
    //if(d->rows == 2)
    //    undecided = norm ? (deltaXY > 0.99 && deltaYX > 0.99) : false;
    if(out.second == 0.0 || undecided){
        output += "\nResult: X--Y with epsilon = " + std::to_string(epsilon);
        direction = "--";
        out.first = 0;
    }else if(deltaXY < deltaYX){
        output += "\nResult: X->Y with epsilon = " + std::to_string(epsilon);
        direction = "->";
        out.first = 1;
    }else{
        output += "\nResult: X<-Y with epsilon = " + std::to_string(epsilon);
        direction = "<-";
        out.first = -1;
    }
    output += "\nRuntime: " +std::to_string(runtime) + " second\n";
    
    if(deltaXY == 1.0 || deltaYX == 1.0){
        output += "CAUTION: X or Y did not compress at all -- maybe the resolution did not suit the data.\n";
    }
    std::cout << "---------------------------\n" << output;
    if(outfolder != NULL && !append && !bagging){
        myfile.open((std::string(outfolder) + "summary.txt").c_str());
        if(myfile.is_open()){
            myfile << "Execution:\t";
            for(int i = 0; i < argc; i++){
                myfile << std::string(argv[i]) + " ";
            }
            myfile << "\n" + output;
        }
        myfile.close();
    }
    if(append && !bagging){
        myfile.open((std::string(outfolder) + "results.tab").c_str(), std::ios_base::app);
        if(myfile.is_open()){
            myfile << std::string(append) + "\t" + std::to_string(epsilon) + "\t" + direction + "\t" + std::to_string(runtime) + "\n";
        }
        myfile.close();
    }
    // delete dot files
    deleteDotFiles(result2);
    deleteDotFiles(result3);
    return out;
}

int main(int argc, char** argv){
    
    // init random numbers
    srand(233314115);
    
    static struct option longopts[] = {
        {"out",             required_argument,  NULL, 'o'},
        {"in",              required_argument,  NULL, 'i'},
        {"type",            required_argument,  NULL, 't'},
        {"precision",       required_argument,  NULL, 'p'},
        {"resolution",      required_argument,  NULL, 'r'},
        {"delimiter",       required_argument,  NULL, 'd'},
        {"minleaveS",       required_argument,  NULL, 'l'},
        {"maxheight",       required_argument,  NULL, 'm'},
        {"score",           required_argument,  NULL, 's'},
        {"x",               required_argument,  NULL, 'x'},
        {"bagging",         required_argument,  NULL, 'b'},
        {"test",            required_argument,  NULL, 'y'},
        {"causal",          no_argument,        NULL, 'c'},
        {"normalize",       no_argument,        NULL, 'n'},
        {"append",          no_argument,        NULL, 'a'},
        {"help",            no_argument,        NULL, 'h'},
        { NULL,             0,                  NULL,  0 }
    };
    
    const char *inname = NULL;
    const char *outfolder = NULL;
    int precision = 4;
    double resolution = 1;
    int xColsSmallerThan = -1;
    int maxHeight = -1;
    bool normalize = false;
    bool causal = true;
    bool ergo = false;
    const char *append = NULL;
    char type = 'X';
    std::string del = "\t";
    double minL = 0.0;
    double bagging = 1;
    int testnr = -1;
    std::set_new_handler(outofmemory);
    bool scoreNormalization = true;
    
    int ch;
    while ((ch = getopt_long(argc, argv, "o:i:b:t:m:p:d:l:y:s:r:x:a:nch", longopts, NULL)) != -1) {
        switch (ch) {
            case 'h':
                printf("Usage: %s -i <input file (tab delimited)> [options]\n", argv[0]);
                printf("  -h    print this help\n");
                printf("  -o    pass output folder in which dot files are outputted (already existing one)\n");
                printf("  -s    pass the scoring method ('1' -- Crack_{\\Delta} (default), '3' -- Crack_{NCI}\n");
                printf("  -t    pass type if all attributes have the same type ('d' .dat-file, 'b' binary, 'c' categorial, 'o' ordinal, 'i' integer, 'n' numeric) -- if not passed it is assumed that the header contains the correct type for each column\n");
                //printf("  -p    precision of the output files (digits behind the floating point -- default: 4)\n");
                printf("  -d    delimiter for text file (default: \\t)\n");
                //printf("  -c    do causal inference -- define cutoff for X and Y with '-x'\n");
                printf("  -x   Defines the index until X goes. Ex. '-x 12' means X includes columns 1 to 12 and Y clumns 13 to M.).\n");
                return 0;
                break;
            case 'i':
                inname = optarg;
                break;
            case 'o':
                outfolder = optarg;
                break;
            case 'l':
                minL = atof(optarg);
                if(minL < 0.0 || minL > 1.0){
                    minL = 0.0;
                    std::cout << "l must be in [0,1]... it was now set to 0.0\n";
                }
                break;
            case 'r':
                resolution = atof(optarg);
                if(resolution <= 0.0 || resolution > 1.0){
                    resolution = 1.0;
                    std::cout << "r must be in (0.0,1.0]... it was now set to 1.0\n";
                }
                break;
            case 'd':
                del = optarg;
                break;
            case 'x':
                xColsSmallerThan = atoi(optarg);
                if(xColsSmallerThan < 1){
                    printf("ERROR: -c argument has to fall into the interval [1, #Attributes - 1] -- the passed argument was '%s'\n", optarg);
                    return 1;
                }
                break;
            case 'c':
                causal = true;
                break;
            case 'n':
                normalize = true;
                break;
            case 'a':
                append = optarg;
                break;
            case 'p':
                precision = atoi(optarg);
                if(precision < 0){
                    precision = 4;
                    std::cout << "Precision (-p) has to be larger than 0... default (4) is used!\n";
                }
                break;
            case 't':
                type = optarg[0];
                if(!(type == 'b' || type == 'c' || type == 'o' || type == 'i' || type == 'n' || type == 'd')){
                    printf("ERROR: -t flag must be choosen from ('d' .dat-file, 'b' binary, 'c' categorial, 'o' ordinal, 'i' integer, 'n' numeric) -- the passed argument was '%s'\n", optarg);
                    return 1;
                }
                break;
            case 'b':
                bagging = atoi(optarg);
                if(bagging < 1){
                    printf("ERROR: -b flag must be passed with the amount of runs. This number is not allowed to be smaller than 1! You passed %s.'\n", optarg);
                    return 1;
                }
                break;
            case 'y':
                testnr = atoi(optarg);
                break;
            case 'm':
                maxHeight = atoi(optarg);
                break;
            case 's':
                int method = atoi(optarg);
                if(method == 2){
                    ergo = false;
                    scoreNormalization = false;
                }else if(method == 3){
                    ergo = true;
                    scoreNormalization = true;
                }else if(method == 4){
                    ergo = true;
                    scoreNormalization = false;
                }else{
                    ergo = false;
                    scoreNormalization = true;
                }
                break;
        }
    }
    
    
    if (inname == NULL && testnr == -1) {
        fprintf(stderr, "Missing input file! Add '-h' flag for help.\n");
        return 1;
    }
    if (outfolder != NULL) {
    	if(!dirExists(outfolder)){
            fprintf(stderr, "Output folder does not exist! Add '-h' flag for help.\n");
            return 1;
        }
    }
    std::cout << "Read data...\n";
    Data *d;
    if(type == 'X')
        d = new Data(inname, del, minL, resolution, xColsSmallerThan, causal, normalize);
    else
        d = new Data(inname, type, del, minL, resolution, xColsSmallerThan, causal, normalize);
    if(maxHeight > 0){
        d->maxHeight = maxHeight;
    }
    std::cout << "Run pack...\n";
    if(xColsSmallerThan == -1){
        pns::Result result = RACK::runGreedyAlgorithm(d, precision, -1);
        if (outfolder != NULL) {
            printOutput(result, outfolder);
        }
    }else{
        if(causal){
            int causalDirection = 0;
            std::map<int,int> counts;
            counts[-1] = 0;
            counts[0] = 0;
            counts[1] = 0;
            std::string cd = "--";
            double t1 = getTime();
            for(int i = 0; i < bagging; i++){
                if(i == 0){
                    if(!ergo)
                        causalDirection = runCracko(d, precision, append, outfolder, bagging > 1, scoreNormalization, argc, argv).first;
                    else
                        causalDirection = runCracker(d, precision, append, outfolder, bagging > 1, scoreNormalization, argc, argv).first;
                }else{
                    Data* dd = new Data(d, randomSample(d->rows, d->rows));
                    if(!ergo)
                        causalDirection = runCracko(dd, precision, append, outfolder, bagging > 1, scoreNormalization, argc, argv).first;
                    else
                        causalDirection = runCracker(dd, precision, append, outfolder, bagging > 1, scoreNormalization, argc, argv).first;
                    counts[causalDirection]++;
                    delete dd;
                }
            }
            double fraction = 0.0;
            if(counts[-1] == counts[1] || (counts[0] > counts[-1] && counts[0] > counts[1])){
                fraction = (double)counts[1] / (double) bagging;
            }else{
                if(counts[-1] > counts[1]){
                    cd = "<-";
                    fraction = (double)counts[-1] / (double) bagging;
                }else{
                    cd = "->";
                    fraction = (double)counts[1] / (double) bagging;
                }
            }
            if(bagging > 1){
                double runtime = getTime() - t1;
                std::string ostring =  std::to_string(fraction) + "\t" + cd + "\t" + std::to_string(runtime) + "\n";
                std::cout << ostring;
                if(append && outfolder != NULL){
                    std::ofstream myfile;
                    myfile.open((std::string(outfolder) + "results.tab").c_str(), std::ios_base::app);
                    if(myfile.is_open()){
                        myfile << std::string(append) + "\t" + ostring;
                    }
                    myfile.close();
                }
            }
        }else{
            pns::Result result = RACK::runGreedyAlgorithm(d, precision, 3);
            if (outfolder != NULL) {
                printOutput(result, outfolder);
            }
        }
    }
    delete d;
}
