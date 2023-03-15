#ifndef TREE_H
#define TREE_H

#include <vector>
#include <utility>
#include <set>

#include "defs.h"
#include "data.h"


class Tree {
    
public:
    Tree(int, std::vector<int>, Data*, int);
    ~Tree();
    void deleteNode(std::vector< pns::node* >);
    
    //std::vector<int> getDependencies();
    double getScore();
    int getPrimeAttribute();
    double decreaseScore(double);
    pns::node* getTree();
    std::vector<int> getCandidates();
    void removeCandidates(std::set<int>);
    bool isCandidate(int);
    std::set<int> dependencies;
    
    pns::dot_file* treeToDot(int);
    bool canSplit;    
    /**
     * Input:
     * 1: contains a list over index sets (list) into which the node was split
     * 2: cut points; list of upper bounds splitting the parent node (size = size of first argument - 1)
     * 3: ref to node that should be split
     * 4: attribute identifier on which the split is based
     *
     * Tasks:
     * assign cut-point list to input node, update its attribute identifier; create children with their
     * corresponding index sets, the old attribute identifier and the given parent node, update their path;
     * add new children to parent node; add new attribute to dependencies
     **/
    void splitNodePack(std::vector< std::vector< int > >, std::vector< double >, pns::node*, int);
    int splitNodeGenerate(std::vector< std::vector< int > >, std::vector< double >, pns::node*, int);
    double getTrivialScore();
    void addOtherOne(int);
protected:
    std::vector<int> candidates;
    pns::node *tree;
    void removeCandidate(int pos);
    Data *data;
    int primeAttribute; // attribute which is encoded by this tree
    double score;
    double trivialScore;
};

#endif
