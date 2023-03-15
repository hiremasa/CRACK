#include <algorithm>
#include <iostream>

#include <map>

#include "forest.h"
#include "score.h"

Forest::Forest(Data *d){
    data = d;
    for(int i = 0; i < d->numA; i++){
        mapping.push_back(-1);  // later to be filled with the actual position of the tree in the forest
        edges.push_back(std::set<int>());
        children.push_back(std::set<int>());
    }
}

Forest::~Forest(){
    while(forest.size() > 0){
        Tree *t = forest.back();
        forest.pop_back();
        delete t;
    }
}

void Forest::addDAG(int p, int c){
    int parent = mapping[p];
    int child = mapping[c];
    if(parent == -1){
        edges[child].insert(p);
    }else{
        edges[child].insert(parent);
        std::set<int> dummy = children[child];
        dummy.insert(child);
        updateDAG(parent, dummy);
    }
}

void Forest::updateDAG(int node, std::set<int> childs){
    // mapping not necessary since only called internally
    if(!pns::isSubset(childs, children[node])){
        children[node] = pns::joinSets(childs, children[node]);
        for (std::set<int>::iterator it = edges[node].begin(); it != edges[node].end(); ++it)
            updateDAG(*it, childs);
    }
}

void Forest::removeCandidates(){
    for(int i = 0; i < forest.size(); i++)
        forest[i]->removeCandidates(children[i]);
}

void Forest::createTrivialTree(int attribute, int whichOne){
    Tree *trivial = new Tree(attribute, data->allIndices, data, whichOne);
    mapping[attribute] = forest.size();
    forest.push_back(trivial);
}

double Forest::getForestScore(){
    return forestScore(data, forest);
}

pns::candidate Forest::findCandidateSplit(int primeAttribute){
    Tree *currentT = forest[mapping[primeAttribute]];
    std::vector< pns::node* > stack;
    stack.push_back(currentT->getTree());
    pns::candidate cand;
    cand.init();
    if(!currentT->canSplit){
        return cand;
    }
    currentT->canSplit = false; // set true as soon as one split was found
    while(stack.size() > 0){
        pns::node *current = stack.back();
        stack.pop_back();
        for(int i = 0; i < current->children.size(); i++){
            stack.push_back(current->children[i]);
        }
        if(current->isLeaf()){
            pns::candidate candL;
            candL.init();
            if(!current->canSplit)
                continue;

            currentT->canSplit = true;
            // check if leaf has already cached costs and if remains DAG
            if(current->splitOnSet && currentT->isCandidate(current->nextSplitOn)){
                // set cached leave
                candL.set(current->change, current, current->nextSplitOn, current->ic);
            }else{
                // search for attribute to split on
                std::vector<int> cands = currentT->getCandidates();
                for(int j = 0; j < cands.size(); j++){
                    int i = cands[j];
                    int pcount = pns::countIndex(current->path, i);
                    if(pcount > 0)
                        continue;
                    if(data->maxHeight != -1 && current->path.size() >= (data->maxHeight))
                        continue;
                    pns::IndCut* cc = splitNode(data, current, i);
                    int count = 0;
                    for(std::vector<int> vec : cc->first){
                        if(vec.size() > 0)
                            count++;
                    }
                    // small hack to accept regression
                    if(cc->second.first[0] == -1.0 && cc->first.size() == 1)
                        count = 2;
                    if(cc->second.second < candL.change && count > 1){
                        candL.set(cc->second.second, current, i, cc);
                    }else{
                        delete cc;
                    }
                }
                if(data->nci){
                    candL.change /= currentT->getTrivialScore();
                }
                // cache costs for node
                if(candL.splitOn != -1){
                    current->splitOnSet = true;
                    current->nextSplitOn = candL.splitOn;
                    current->change = candL.change;
                    current->ic = candL.ic;
                }else{
                    current->canSplit = false;
                }
            }
            // update overall best split on tree
            if(candL.change < cand.change){
                cand.set(candL);
            }
        }
    }
    return cand;
}

void Forest::split(pns::candidate cand){
    Tree *actual = forest[mapping[cand.toSplit->attribute]];
    actual->splitNodePack(cand.ic->first, cand.ic->second.first, cand.toSplit, cand.splitOn);
}

std::vector< pns::dot_file* > Forest::forestToDot(int precision){
    std::vector < pns::dot_file* > outp;
    for(int i = 0; i < forest.size(); i++){
        outp.push_back(forest[i]->treeToDot(precision));
    }
    return outp;
}

void Forest::addOtherOnes(int whichOne){
    for(int i = 0; i < forest.size(); i++){
        forest[i]->addOtherOne(whichOne);
    }
}

int Forest::size(){
    return forest.size();
}

void Forest::printStats(){
    std::map<int,int> stats;
    for(int i = 0; i < forest.size(); i++){
        Tree *currentT = forest[i];
        std::vector< pns::node* > stack;
        stack.push_back(currentT->getTree());
        pns::candidate cand;
        cand.init();
        while(stack.size() > 0){
            pns::node *current = stack.back();
            stack.pop_back();
            for(int i = 0; i < current->children.size(); i++){
                stack.push_back(current->children[i]);
            }
            if(current->isLeaf()){
                stats[current->indices.size()]++;
            }
        }
    }
    std::cout << "Leaf sizes:\n";
    int maxc = (int) round((double) data->rows * 0.01);
    int sumc = 0;
    int sumt = 0;
    for(std::map<int,int>::iterator it = stats.begin(); it != stats.end(); ++it){
        if(it->first < maxc){
            sumc += it->second;
        }
        sumt += it->second;
    }
    std::cout << (double)sumc / (double)sumt << std::endl;
}
