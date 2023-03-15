#ifndef FOREST_H
#define FOREST_H

#include <vector>
#include <utility>
#include <set>

#include "tree.h"
#include "data.h"
#include "defs.h"

class Forest {
    
public:
    Forest(Data*);
    ~Forest();
    
    /**
     * 1: Attribute the split is based on
     * 2: Attribute the tree belongs to (prime attribute)
     **/
    //bool remainsDAG(int,int);
    void createTrivialTree(int,int);
    pns::candidate findCandidateSplit(int);
    void split(pns::candidate);
    double getForestScore();
    void addDAG(int,int);
    void removeCandidates();
    // convention: first line is name of file -- pass precision
    std::vector< pns::dot_file* > forestToDot(int);
    void printStats();
    int size();
    void addOtherOnes(int);
protected:
    std::vector<Tree*> forest;
    std::vector< std::set<int> > edges;
    std::vector< std::set<int> > children;
    std::vector<int> mapping; // maps position of tree to attribute; allowing constant access
    Data *data;
    void updateDAG(int, std::set<int>);
};

#endif
