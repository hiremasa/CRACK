#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <limits>
#include <cmath>

#include "data.h"

Data::Data(const char* file, char type, std::string del, double pp, double res, int xCST, bool c, bool n){
    nci = false;
    causal = c;
    xColsSmallerThan = xCST;
    if(type == 'd'){
        read_dat(file);
    
    }else{
        typeIndicator.push_back(type);
        read_data(file, del);
    }
    std::vector<double> dummyR(numA, res);
    resolution = dummyR;
    minRes = res;
    if(n)
        normalizeData();
    else
        preprocessData();
    setStandards(pp);
}

Data::Data(const char* file, std::string del, double pp, double res, int xCST, bool c, bool n){
    causal = c;
    nci = false;
    xColsSmallerThan = xCST;
    read_data(file, del);
    std::vector<double> dummyR(numA, res);
    resolution = dummyR;
    minRes = res;
    if(n)
        normalizeData();
    else
        preprocessData();
    setStandards(pp);
}

Data::Data(Data* d, std::vector<int> indices){
    numA = d->numA;
    nci = d->nci;
    maxHeight = d->maxHeight;
    rows = indices.size();
    resolution = d->resolution;
    minRes = d->minRes;
    worstCaseDomain = d->worstCaseDomain;
    causal = d->causal;
    minLeafSize = d->minLeafSize;
    binBound = d->binBound;
    xColsSmallerThan = d->xColsSmallerThan;
    cache = d->cache;
    log2faccache = d->log2faccache;
    domain = d->domain;
    typeIndicator = d->typeIndicator;
    internalCall = false;
    for(int i = 0; i < indices.size(); i++){
        allIndices.push_back(i);
        data.push_back(d->data[indices[i]]);
    }
}

Data::Data(int m, int k, int l){
    allIndices = std::vector<int>();
    for(int i = 0; i < m; i++){
        allIndices.push_back(i);
        data.push_back(std::vector<double>((k+l), 0.0));
    }
    numA = (k+l);
    rows = m;
    nci = false;
    std::vector<double> dummyR(numA, 1.0);
    resolution = dummyR;
    minRes = 1.0;
    xColsSmallerThan = k;
    causal = true;
    std::vector<double> dummy(rows, -1.0);
    std::vector<double> dummy2(rows, -1.0);
    cache = dummy;
    log2faccache = dummy2;
    internalCall = false;
    typeIndicator = std::vector<char>();
    std::vector<int> dummyd(numA, 1);
    domain = dummyd;
}

void Data::setStandards(double pp){
    maxHeight = -1;
    std::vector<int> dummy(numA, 0);
    domain = dummy;
    minLeafSize = (int) floor((double) rows * pp);
    if(numA <= xColsSmallerThan){
        std::cout << "The argument passed with -c has to be smaller than #Attributes (" << numA << ")! Currently it is: " << xColsSmallerThan << std::endl;
        exit(1);
    }
    binBound = 0.6;// hardcoded needed for histograms
    double totalMax = -std::numeric_limits<double>::max();
    double totalMin = std::numeric_limits<double>::max();
    double minDomain = totalMin;
    // set resolutions
    for(int j = 0; j < numA; j++){
        char t = getAttributeType(j);
        if(!(t == 'b' || t == 'c' || t == 'o')){
            resolution[j] = getMinDiff(j);
        }
    }
    minRes = totalMin;
    for(int att = 0; att < numA; att++){
        char t = getAttributeType(att);
        if(t == 'b' || t == 'c' || t == 'o'){
            pns::ad* dummyad = getAttributeDistribution(att);
            domain[att]= dummyad->size();
            delete dummyad;
        }else{
            double max = -std::numeric_limits<double>::max();
            double min = std::numeric_limits<double>::max();
            for(int i = 0; i < rows; i++){
                double val = data[i][att];
                if(val < min)
                    min = val;
                if(val > max)
                    max = val;
            }
            domain[att] = ((max - min + 1) * (1.0 / resolution[att]));
            if(max > totalMax)
                totalMax = max;
            if(min < totalMin)
                totalMin = min;
            if(domain[att] < minDomain)
                minDomain = domain[att];
            if(resolution[att] < minRes)
                minRes = resolution[att];
        }
    }
    worstCaseDomain = ((totalMax - totalMin + 1) * (1.0 / minRes));
    fprintf(stdout, "Resolution\tDomain\n");
    for(int j = 0; j < numA; j++){
        fprintf(stdout, "%f\t%d\n", resolution[j], domain[j]);
    }
}

void Data::read_dat(const char* filename){
    internalCall = false;
    
    std::cout << "Extracting data matrix from: " << filename << "\n";
    
    std::ifstream f;
    std::string line;
    rows = 0;
    
    f.open(filename);
    if (f.is_open()) {
        while (std::getline(f, line)) {
            if(line.length() > 0){
                std::vector<std::string> line_elems = splitLine(line, " ");
                data.push_back(strToDouble(line_elems));
                rows++;
            }
        }
    }
    // find min and max
    int max = 0;
    int min = std::numeric_limits<int>::max();
    for(int i = 0; i < data.size(); i++){
        for(int j = 0; j < data[i].size(); j++){
            if(data[i][j] > max)
                max = data[i][j];
            if(data[i][j] < min)
                min = data[i][j];
        }
    }
    numA = (max + 1) - min;
    for(int i = 0; i < data.size(); i++){
        std::vector<double> dummy(numA, 0.0);
        for(int j = 0; j < data[i].size(); j++){
            dummy[data[i][j] - min] = 1;
        }
        data[i] = dummy;
    }
    // set type
    for(int i = 0; i < numA; i++){
        typeIndicator.push_back('b');
    }
    for(int i = 0; i < rows; i++){
        allIndices.push_back(i);
    }
    std::vector<double> dummy(rows, -1.0);
    std::vector<double> dummy2(rows, -1.0);
    cache = dummy;
    log2faccache = dummy2;
    std::cout << numA << " attributes and " << rows << " rows read.\n";
}

void Data::read_data(const char* filename, std::string delimiter){
    internalCall = false;
    
    std::cout << "Extracting data matrix from: " << filename << "\n";
    
    std::ifstream f;
    std::string line;
    
    rows = 0;
    bool storeType = typeIndicator.size() == 0;
    bool firstRow = true;
    f.open(filename);
    if (f.is_open()) {
        while (std::getline(f, line)) {
            if(line.length() > 0){
                std::vector<std::string> line_elems = splitLine(line, delimiter);
                if(delimiter == " " && line_elems.size() <= 1){
                    line_elems = splitLine(line, "\t");
                }
                if(firstRow){
                    firstRow = false;
                    numA = line_elems.size();
                    if(storeType){
                        // type in first row
                        typeIndicator = strToChar(line_elems);
                    }else{
                        // duplicate known type to all rows
                        for(int i = 1; i < numA; i++){
                            typeIndicator.push_back(typeIndicator[0]);
                        }
                        // read first line
                        data.push_back(strToDouble(line_elems));
                        rows++;
                    }
                }else{
                    // read current line
                    assert(line_elems.size() == numA);
                    data.push_back(strToDouble(line_elems));
                    rows++;
                }
            }
        }
    }
    for(int i = 0; i < rows; i++){
        allIndices.push_back(i);
    }
    std::vector<double> dummy(rows, -1.0);
    std::vector<double> dummy2(rows, -1.0);
    cache = dummy;
    log2faccache = dummy2;
    std::cout << numA << " attributes and " << rows << " rows read.\n";
}

double Data::getMinDiff(int attr){
    std::vector<double> values;
    for(int i = 0; i < rows; i++){
        values.push_back(data[i][attr]);
    }
    std::sort(values.begin(), values.end());
    double diff = 1.0;
    std::vector<double> diffs;
    for(int i = 1; i < rows; i++){
        double curr_diff = values[i] - values[i-1];
        if(!pns::isZero(curr_diff)){
            diffs.push_back(curr_diff);
        }
    }
    std::sort(diffs.begin(), diffs.end());
    if(diffs.size() > 0){
        int pos = std::floor((double) diffs.size() * 0.1);
        diff = diffs[pos];
    }
    return diff;
}

void Data::normalizeData(){
    for(int j = 0; j < numA; j++){
        if(getAttributeType(j) == 'b' || getAttributeType(j) == 'c')
            continue;
        double max = -std::numeric_limits<double>::max();
        double min = std::numeric_limits<double>::max();
        for(int i = 0; i < rows; i++){
            double v = data[i][j];
            if(v < min)
                min = v;
            if(v > max)
                max = v;
        }
        for(int i = 0; i < rows; i++){
            double v = data[i][j];
            data[i][j] = ((v - min) / (max - min));
        }
    }
}

void Data::preprocessData(){
    for(int j = 0; j < numA; j++){
        if(getAttributeType(j) == 'b' || getAttributeType(j) == 'c')
            continue;
        resolution[j] = getMinDiff(j);
    }
}

pns::ad* Data::getAttributeDistribution(int attribute){
    internalCall = true;
    pns::ad* dummy = getRestrictedAttributeDistribution(attribute, std::vector<int>());
    internalCall = false;
    return dummy;
}

bool Data::allowedToRegress(){
    return !(causal && numA <= 2);
}

double Data::getTransformedValue(pns::node* n, int pos){
    double v = data[pos][n->attribute];
    double vold = v;
    v -= n->alpha;
    for(int j = 0; j < n->beta.size(); j++){
        double x = data[pos][n->beta_ref[j]];
        v -= n->beta[j] * x;
        v -= n->gamma[j] * (x * x);
    }
    return v;
}

double Data::getLog2Fac(int n){
    if(log2faccache[n] == -1){
        log2faccache[n] = pns::log2fac(n);
    }
    return log2faccache[n];
}

pns::ad* Data::getRestrictedAttributeDistribution(int attribute, std::vector<int> indices){
    if(internalCall){
        indices = allIndices;
    }
    pns::ad* toReturn = new pns::ad();
    for(int i = 0; i < indices.size(); i++){
        int pos = indices[i];
        (*toReturn)[data[pos][attribute]] += 1;
    }
    return toReturn;
}

std::vector<std::string> Data::splitLine(std::string line, std::string del){
    std::vector<std::string> elements;
    size_t pos = 0;
    while ((pos = line.find(del)) != std::string::npos) {
        std::string dummy = line.substr(0, pos);
        if(dummy != " " && dummy != "\t" && dummy != "")
            elements.push_back(dummy);
        line.erase(0, pos + del.length());
    }
    if(line != " " && line != "\t" && line != "")
        elements.push_back(line);
    return elements;
}

std::vector<double> Data::strToDouble(std::vector<std::string> elements){
    std::vector<double> row;
    for(int i = 0; i < elements.size(); i++){
        double dummy = atof(elements[i].c_str());
        row.push_back(dummy);
    }
    return row;
}

std::vector<char> Data::strToChar(std::vector<std::string> elements){
    std::vector<char> row;
    for(int i = 0; i < elements.size(); i++){
        if(elements[i].length() != 1)
            fprintf(stdout, "Error while reading text file. Please check if formatting is consistent and if you passed or marked the type of the columns correctly!\n");
        assert(elements[i].length() == 1);
        row.push_back(elements[i][0]);
    }
    return row;
}

char Data::getAttributeType(int att){
    return typeIndicator[att];
}
