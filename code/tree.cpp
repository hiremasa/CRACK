#include <iostream>
#include <assert.h>
#include <cmath>
#include <limits>

#include "tree.h"
#include "score.h"

/**
 * whichOne:    0 candidates are X (first part)
 *              1 candidates are Y (second part)
 *              2 candidates are both (Y|X) or vice versa, or everything
 **/
Tree::Tree(int attribute, std::vector<int> indices, Data *d, int whichOne){
    data = d;
    primeAttribute = attribute;
    tree = new pns::node();
    tree->attribute = primeAttribute;
    tree->indices = indices;
    candidates = std::vector<int>();
    pns::ad* dummyad = d->getAttributeDistribution(attribute);
    tree->categories = dummyad->size();
    delete dummyad;
    tree->canSplit = tree->categories > 1; // O(m)
    if(tree->canSplit){
        int start = whichOne == 1 ? data->xColsSmallerThan : 0;
        int end = whichOne == 0 ? data->xColsSmallerThan : data->numA;
        for(int i = start; i < end; i++){
            if(i != attribute)
                candidates.push_back(i);
        }
    }
    canSplit = tree->canSplit;
    tree->hasLeafData = false;
    tree->parent = NULL;
    tree->alpha = 0.0;
    tree->beta = std::vector<double>();
    tree->gamma = std::vector<double>();
    tree->beta_ref = std::vector<int>();
    tree->isRegressor = false;
    tree->multiwaySplit = false;
    tree->duplicates = false;
    trivialScore = treeScore(d, tree, primeAttribute);
    score = trivialScore;
    dependencies = std::set<int>();
}

void Tree::addOtherOne(int whichOne){
    int start = whichOne == 0 ? data->xColsSmallerThan : 0;
    int end = whichOne == 1 ? data->xColsSmallerThan : data->numA;
    for(int i = start; i < end; i++){
        if(i != primeAttribute)
            candidates.push_back(i);
    }
    tree->canSplit = true;
    std::vector< pns::node* > stack;
    stack.push_back(tree);
    while(stack.size() > 0){
        pns::node *current = stack.back();
        stack.pop_back();
        for(int i = 0; i < current->children.size(); i++){
            stack.push_back(current->children[i]);
        }
        if(current->isLeaf()){
            current->canSplit = true;
            current->splitOnSet = false;
            current->nextSplitOn = 0;
        }
    }
}

Tree::~Tree(){
    std::vector<pns::node*> stack;
    stack.push_back(tree);
    deleteNode(stack);
}

void Tree::deleteNode(std::vector<pns::node*> stack){
    while(stack.size() > 0){
        pns::node *n = stack.back();
        stack.pop_back();
        deleteNode(n->children);
        delete n->ic;
        delete n;
    }
}

bool Tree::isCandidate(int cand){
    return pns::containsIndex(candidates, cand) != -1;
}

void Tree::removeCandidates(std::set<int> toRemove){
    for (std::set<int>::iterator it = toRemove.begin(); it != toRemove.end(); ++it){
        int pos = pns::containsIndex(candidates, *it);
        if(pos != -1)
            removeCandidate(pos);
    }
}

void Tree::removeCandidate(int pos){
    candidates.erase(candidates.begin() + pos);
}

std::vector<int> Tree::getCandidates(){
    return candidates;
}

double Tree::getScore(){
    return score;
}

int Tree::getPrimeAttribute(){
    return primeAttribute;
}

pns::node* Tree::getTree(){
    return tree;
}

double Tree::decreaseScore(double change){
    score += change;
    return score;
}

double Tree::getTrivialScore(){
    return trivialScore;
}

void Tree::splitNodePack(std::vector< std::vector< int > > indices, std::vector< double > cps, pns::node* toSplit, int splitOnID){
    assert(cps.size() > 0);
    // buffer id of node to split
    int oldID = toSplit->attribute;
    toSplit->attribute = splitOnID;
    toSplit->cut_points = cps;
    toSplit->multiwaySplit = false;
    if(cps[0] != -1.0 || cps.size() < 3){
        toSplit->isRegressor = false;
        if(indices.size() > 2){
            toSplit->multiwaySplit = true;
        }else{
            assert(indices.size() == (cps.size() + 1));
        }
    }else{
        assert(indices.size() == 1 && cps.size() >= 3);
        toSplit->isRegressor = true;
    }
    dependencies.insert(splitOnID);
    std::vector<int> p = toSplit->path;
    p.push_back(splitOnID);
    for(int i = 0; i < indices.size(); i++){
        if(indices[i].size() < 1)
            continue;
        pns::node *current = new pns::node();
        current->attribute = oldID;
        current->parent = toSplit;
        current->indices = indices[i];
        current->path = p;
        current->hasLeafData = false;
        current->leafData = 0.0;
        pns::ad* dummyad = data->getRestrictedAttributeDistribution(current->attribute, current->indices);
        current->categories = dummyad->size();
        delete dummyad;
        current->canSplit = tree->categories > 1;
        current->canSplit = current->canSplit && current->indices.size() >= data->minLeafSize;
        current->splitOnSet = false;
        toSplit->children.push_back(current);
        current->alpha = toSplit->alpha;
        current->beta = toSplit->beta;
        current->beta_ref = toSplit->beta_ref;
        current->gamma = toSplit->gamma;
        current->isRegressor = false;
        current->multiwaySplit = false;
        current->duplicates = toSplit->multiwaySplit && i != 0;
        if(cps[0] == -1.0 && !toSplit->multiwaySplit){
            current->alpha += cps[1];
            current->beta.push_back(cps[2]);
            current->beta_ref.push_back(splitOnID);
            if(cps.size() == 4){
                current->gamma.push_back(cps[3]);
            }else{
                current->gamma.push_back(0.0);
            }
            rvLeafData(data, current);
        }
    }
}

int Tree::splitNodeGenerate(std::vector< std::vector< int > > indices, std::vector< double > cps, pns::node* toSplit, int splitOnID){
    assert(cps.size() > 0);
    // buffer id of node to split
    int oldID = toSplit->attribute;
    toSplit->attribute = splitOnID;
    toSplit->cut_points = cps;
    if(cps[0] != -1.0 || cps.size() < 3){
        toSplit->isRegressor = false;
        assert(indices.size() == (cps.size() + 1));
    }else{
        assert(indices.size() == 1 && cps.size() >= 3);
        toSplit->isRegressor = true;
    }
    std::vector<int> p = toSplit->path;
    p.push_back(splitOnID);
    int count = -1;
    for(int i = 0; i < indices.size(); i++){
        if(indices[i].size() == 0)
            continue;
        pns::node *current = new pns::node();
        current->attribute = oldID;
        current->parent = toSplit;
        current->indices = indices[i];
        current->path = p;
        current->hasLeafData = false;
        current->leafData = 0.0;
        pns::ad* dummyad = data->getRestrictedAttributeDistribution(current->attribute, current->indices);
        current->categories = dummyad->size();
        delete dummyad;
        current->canSplit = tree->categories > 1;
        current->canSplit = current->canSplit && current->indices.size() >= data->minLeafSize;
        current->splitOnSet = false;
        toSplit->children.push_back(current);
        current->alpha = toSplit->alpha;
        current->beta = toSplit->beta;
        current->beta_ref = toSplit->beta_ref;
        current->gamma = toSplit->gamma;
        current->isRegressor = false;
        current->multiwaySplit = false;
        if(cps[0] == -1.0){
            current->alpha += cps[1];
            current->beta.push_back(cps[2]);
            current->beta_ref.push_back(splitOnID);
            if(cps.size() == 4){
                current->gamma.push_back(cps[3]);
            }else{
                current->gamma.push_back(0.0);
            }
            rvLeafData(data, current);
        }
        count++;
    }
    return count;
}

pns::dot_file* Tree::treeToDot(int precision){
    pns::dot_file* file = new pns::dot_file();
    std::string id = std::to_string(primeAttribute + 1) + "S";
    // insert file name
    file->push_back("subgraph cluster" + id + " {\n");
    file->push_back("\tlabel = \"" + std::to_string(primeAttribute + 1) + " (" + data->getAttributeType(primeAttribute) + "):\";\n");
    int index = 0;
    std::vector< pns::node* > stack;
    std::vector< int > pos_stack;
    pos_stack.push_back(0);
    stack.push_back(tree);
    while(stack.size() > 0){
        pns::node *current = stack.back();
        stack.pop_back();
        int pos = pos_stack.back();
        pos_stack.pop_back();
        if(current->isLeaf()){
            // check distribution
            pns::ad* mapping = data->getRestrictedAttributeDistribution(current->attribute, current->indices);
            file->push_back("\ta" + id + std::to_string(pos) + " [shape=box, label=\"#" + std::to_string(current->indices.size()) + "\\n");
            int cp_pos = 0;
            int count = 0;
            int total = current->indices.size();
            char type = data->getAttributeType(current->attribute);
            if(type == 'c' || type == 'b'){
                for(pns::ad::iterator it = mapping->begin(); it != mapping->end(); it++) {
                    double key = it->first;
                    double count = (double) it->second;
                    // round to precision digits
                    double prob = count / (double) total;
                    file->push_back("p(a" + std::to_string((current->attribute + 1)) + " = " + pns::pruneDouble(key,precision) + ") = " + pns::pruneDouble(prob, precision) + "\\n");
                }
            }else{
                double t = 0.0;
                double tt = 0.0;
                double N = (double) current->indices.size();
                for(int i = 0; i < current->indices.size(); i++){
                    int posp = current->indices[i];
                    double v = data->getTransformedValue(current, posp);
                    t += v;
                    tt += (v*v);
                }
                double mu = (double)t / N;
                double sse = tt - (((double)1/N) * (t*t));
                double var = sqrt(sse / N);
                if(var != var)
                    var = 0.0;
                if(current->beta.size() > 0){
                    std::string line = pns::pruneDouble(current->alpha + mu, precision);
                    int variables = 0;
                    for(int q = 0; q < current->beta.size(); q++){
                        if(variables % 2 == 1)
                            line += "\\n";
                        line += " + " + pns::pruneDouble(current->beta[q], precision) + "a" + std::to_string(current->beta_ref[q] + 1);
                        variables++;
                        if(current->gamma[q] != 0.0){
                            if(variables % 2 == 1)
                                line += "\\n";
                            line += " + " + pns::pruneDouble(current->gamma[q], precision) + "a" + std::to_string(current->beta_ref[q] + 1) + "&#xb2;";
                            variables++;
                        }
                    }
                    file->push_back(line + "\\n");
                }else{
                    file->push_back("&#956; = " + pns::pruneDouble(mu,precision) + "\\n");
                }
                file->push_back("&#963; = " + pns::pruneDouble(var,precision) + "\\n");
            }
            file->push_back("\"];\n");
            delete mapping;
        }else{
            // create nodes, indices and arcs
            file->push_back("\ta" + id + std::to_string(pos) + " [label=\"a" + std::to_string((current->attribute + 1)) + "\", shape=");
            if(current->isRegressor)
                file->push_back("parallelogram];\n");
            else
                file->push_back("diamond];\n");
            for(int i = 0; i < current->children.size(); i++){
                index++;
                pos_stack.push_back(index);
                pns::node *child = current->children[i];
                stack.push_back(child);
                if(current->isRegressor){
                    file->push_back("\ta" + id + std::to_string(pos) + " -> a" + id + std::to_string(index) + " [label=\"&#945; = " + pns::pruneDouble(current->cut_points[1], precision) + "\"];\n");
                }else{
                    char t = data->getAttributeType(current->attribute);
                    if(i < current->cut_points.size()){
                        if(t == 'b'){
                            file->push_back("\ta" + id + std::to_string(pos) + " -> a" + id + std::to_string(index) + " [label=\"= 0\"];\n");
                        }else if(t == 'c'){
                            file->push_back("\ta" + id + std::to_string(pos) + " -> a" + id + std::to_string(index) + " [label=\"&#8800; " + pns::pruneDouble(current->cut_points[i], precision) + "\"];\n");
                        }else{
                        file->push_back("\ta" + id + std::to_string(pos) + " -> a" + id + std::to_string(index) +   " [label=\"< " + pns::pruneDouble(current->cut_points[i], precision) + "\"];\n");
                        }
                    }else{
                        if(t == 'b'){
                            file->push_back("\ta" + id + std::to_string(pos) + " -> a" + id + std::to_string(index) +    " [label=\"= 1\"];\n");
                        }else if(t == 'c'){
                            file->push_back("\ta" + id + std::to_string(pos) + " -> a" + id + std::to_string(index) + " [label=\"= " + pns::pruneDouble(current->cut_points[(i-1)], precision) + "\"];\n");
                        }else{
                            file->push_back("\ta" + id + std::to_string(pos) + " -> a" + id + std::to_string(index) + " [label=\"&#8805; " + pns::pruneDouble(current->cut_points[(i-1)],precision) + "\"];\n");
                        }
                    }
                }
            }
        }
    }
    file->push_back("\tgraph[style=dotted];\n");
    file->push_back("}\n");
    return file;
}
