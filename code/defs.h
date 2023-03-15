#ifndef DEFS_H
#define DEFS_H

#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <set>

namespace pns { //pack name space
    typedef std::map<double,int> ad; // attribute distribution
    typedef std::vector< std::vector<double> > data_map;
    //IndCut:   first: list of indices
    //          second.first: list of cut-points (if first is -1.0, 2nd is alpha and 3rd is beta)
    //          second.second: score change
    typedef std::pair<std::vector<std::vector<int> >, std::pair< std::vector< double >, double > > IndCut;
    typedef std::vector< std::string > dot_file;
    
    struct node {   // always use new
        int attribute;
        node *parent;
        std::vector< node* > children;
        std::vector< double > cut_points;  // cut points (most probably only one -- for binary split)
        std::vector< int > indices; // indices of the rows covered by this node
        std::vector< int > path;
        double alpha;
        std::vector<double> beta;
        std::vector<double> gamma;
        std::vector<int> beta_ref;
        bool splitOnSet;    // init with false
        int nextSplitOn;    // attribute on which to split
        bool hasLeafData;     // leaf data costs set?
        double leafData;   // leaf data costs
        double change; // the more negative the better
        bool isRegressor;
        bool multiwaySplit;
        bool duplicates;
        bool canSplit;
        int categories;
        IndCut* ic;
        
        bool isLeaf(){
            return children.size() == 0;
        }
    };
    
    struct candidate {
        double change;
        int splitOn;
        node* toSplit;
        IndCut* ic;
        void init(){ change = 0.0; splitOn = -1; toSplit = 0; }
        void set(candidate b){
            change = b.change;
            toSplit = b.toSplit;
            splitOn = b.splitOn;
            ic = b.ic;
        }
        void set(double c, node* t, int s, IndCut* i){
            change = c;
            toSplit = t;
            splitOn = s;
            ic = i;
        }
        void deleteCut(){
            delete ic;
        }
    };
    
    int inline containsIndex(std::vector<int> v, int i){
        int pos = std::find(v.begin(), v.end(), i) - v.begin();
        return pos < v.size() ? pos : -1;
    }
    
    int inline countIndex(std::vector<int> v, int i){
        int count = 0;
        int pos = std::find(v.begin(), v.end(), i) - v.begin();
        while(pos < v.size()){
            count++;
            pos = std::find(v.begin() + pos, v.end(), i) - v.begin();
            pos++;
        }
        return count;
    }
    
    std::string inline pruneDouble(double d, int precision){
        int actual_precision = 0;
        double dummy = d;
        for(int i = 0; i <= precision; i++){
            double intpart;
            actual_precision = i;
            if(modf(dummy, &intpart) == 0.0){
                break;
            }else{
                dummy *= 10;
            }
        }
        std::stringstream stream;
        stream << std::fixed << std::setprecision(actual_precision) << d;
        return stream.str();
    }
    
    void inline printSet(int ID, std::set<int> s){
        std::cout << ID << ": <";
        int elem = 0;
        for (std::set<int>::iterator it = s.begin(); it != s.end(); ++it){
            if(elem == 0)
                std::cout << *it;
            else
                std::cout << ", " << *it;
            elem++;
        }
        std::cout << ">\n";
    }
    
    bool inline isSubset(std::set<int> a, std::set<int> b){
        std::set<int>::iterator it;
        for (it = a.begin(); it != a.end(); ++it){
            if(!(b.find(*it) != b.end()))
                return false;
        }
        return true;
    }
    
    std::set<int> inline joinSets(std::set<int> a, std::set<int> b){
        std::set<int>::iterator it;
        for (it = a.begin(); it != a.end(); ++it)
            b.insert(*it);
        return b;
    }
    
    struct IndValPair {
        int pos;
        int end; // this index is excluded
        double value;
    };
    
    struct Regressor {
        double alpha;
        double beta;
        bool hasGamma;
        double gamma;
        double sse;
    };
    
    struct by_value {
        bool operator()(IndValPair const &a, IndValPair const &b) {
            return a.value < b.value;
        }
    };
    
    double inline log2fac(int n){
        double sum = 0;
        for(int i = 2; i <= n; i++){
            sum += log2(i);
        }
        return sum;
    }
    
    struct Result {
        std::vector< pns::dot_file* > files;
        double initialScore;
        double finalScore;
        double runtime;
        
        double getCompression(){
            return finalScore / initialScore;
        }
        void set(std::vector< pns::dot_file* > f, double iS, double fS, double r){
            files = f;
            initialScore = iS;
            finalScore = fS;
            runtime = r;
        }
    };
    
    double inline log2nChoosek(int n, int k){
        if(k > n || k == 0){
            return 0;
        }else{
            return log2fac(n) - log2fac(k) - log2fac(n-k);
        }
    }
    
    double inline min(int a, int b){
        return a < b ? a : b;
    }
    
    double inline min(double a, double b){
        return a < b ? a : b;
    }
    
    double inline log2N(double z){
        if(z < 1){
            std::cout << "ERROR: log2N not defined for z < 1! --> returned 0\n";
            return 0;
        }else{
            double logstar = log2(z);
            double sum = logstar;
            while(true){
                logstar = log2(logstar);
                if(logstar <= 0){
                    break;
                }else{
                    sum += logstar;
                }
            }
            return sum + log2(2.865064);
        }
    }
    
    double inline absD(double d){
        return d < 0 ? -d : d;
    }
    
    bool inline isZero(double x){
        x = absD(x);
        return x < 0.0000001;
    }
}

#endif
