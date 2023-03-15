#ifndef DATA_H
#define DATA_H

#include <utility>
#include <string>
#include <vector>

#include "defs.h"



class Data {
    
public:
    bool nci;
    int numA; // number of attributes
    int rows; // number of rows per attribute
    std::vector<double> resolution;
    double minRes;
    double worstCaseDomain;
    bool causal;
    int minLeafSize;
    double binBound;
    int xColsSmallerThan;
    std::vector<int> allIndices;
    pns::data_map data;                 // vector over attributes; each attribute contains its corresponding values
    std::vector<double> cache;
    std::vector<double> log2faccache;
    std::vector<int> domain; // set while tree initialization (forest building)
    std::vector<char> typeIndicator;    // 'b' binary, 'c' categorial, 'n' nominal, 'i' integer, 'n' numeric
    
    Data(const char*,char, std::string,double,double,int,bool, bool); // filename, data type, delimiter, in [0,1] percentage of rows to not split anymore
    Data(const char*,std::string,double,double, int,bool, bool);       // filename, delimiter,  first row must consist of the type indicators for mixed types
    Data(Data*, std::vector<int>); // copy bagging sample
    Data(int, int, int);
    pns::ad* getAttributeDistribution(int);  // get complete attribute distribution
    pns::ad* getRestrictedAttributeDistribution(int, std::vector<int>); // get attribute distribution over subspace
    char getAttributeType(int);
    double getLog2Fac(int);
    double getTransformedValue(pns::node*, int);
    bool allowedToRegress();
    void normalizeData();
    void setStandards(double);
    void preprocessData();
    
    // max height is only used for testing
    int maxHeight;
    
protected:
    void read_data(const char*,std::string);
    void read_dat(const char*);
    double getMinDiff(int);
    std::vector<std::string> splitLine(std::string, std::string);
    std::vector<double> strToDouble(std::vector<std::string>);
    std::vector<char> strToChar(std::vector<std::string>);
    
private:
    bool internalCall;
};

#endif
