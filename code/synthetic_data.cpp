#include "synthetic_data.h"

#include <stdlib.h> // rand

void generateRandomData(int k, int l, char type, int domainShift, double splitProbability, double dependencyProbability, Data* d){
    bool sameDomain = domainShift != 0;
    int sprob = round(splitProbability * 100.0);
    int dprob = round(dependencyProbability * 100.0);
    for(int i = 0; i < (k+l); i++){
        // save type of attribute i
        if(type == 'a' || type == 'z'){ //all
            int diffTypes = 3;
            if(type == 'z')
                diffTypes = 2;
            int random = rand() % diffTypes;
            if(random == 0){
                d->typeIndicator.push_back('b');
            }else if(random == 1){
                d->typeIndicator.push_back('c');
            }else{
                d->typeIndicator.push_back('i');
            }
        }else{
            d->typeIndicator.push_back(type);
        }
        std::cout << d->getAttributeType(i) << " ";
        Tree *current = new Tree(i, d->allIndices, d, -1);
        int leaves = 1;
        int start = 0;
        if(i >= k){
            int rdep = rand() % 100 + 1;
            if(rdep > dprob)
                start = k;
        }
        for(int j = start; j < i; j++){ //perform splits
            int r = rand() % 100 + 1;
            if(r <= sprob){ // do split or similar
                leaves += splitLeaf(current, d, i , j, leaves, domainShift);
            }// no split
        }
        if(i < k)
            generateDataFromLeaves(current->getTree(), d, i, 0);
        else
            generateDataFromLeaves(current->getTree(), d, i, domainShift);
        delete current;
    }
    std::cout << std::endl;
}

//Box muller method
double rand_normal(double mean, double stddev){
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;
            
            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

/**
 * returns 1 for split and 0 for regression
 **/
int splitLeaf(Tree* T, Data* d, int a, int c, int children, int domainShift){
    char tA = d->getAttributeType(a);
    char tC = d->getAttributeType(c);
    char type = d->getAttributeType(a);
    int whichChild = rand() % children;
    int domain = 5;
    int sqrdom = round(10 * sqrt(domain));
    int currChild = 0;
    std::vector<pns::node*> stack;
    pns::node *t = T->getTree();
    stack.push_back(t);
    while(stack.size() > 0){
        pns::node *current = stack.back();
        stack.pop_back();
        for(int i = 0; i < current->children.size(); i++){
            stack.push_back(current->children[i]);
        }
        if(current->isLeaf()){
            if(currChild == whichChild){
                std::vector< std::vector< int > > indices;
                std::vector< double > cps;
                pns::node* toSplit = current;
                int splitOnID = c;
                if(tC == 'b'){ // binary
                    std::vector<int> i0;
                    std::vector<int> i1;
                    for(int i = 0; i < current->indices.size(); i++){
                        int pos = current->indices[i];
                        if(d->data[pos][c] == 0)
                            i0.push_back(pos);
                        else
                            i1.push_back(pos);
                    }
                    cps.push_back(1);
                    indices.push_back(i0);
                    indices.push_back(i1);
                    return T->splitNodeGenerate(indices, cps, toSplit, splitOnID);
                }else if(tC == 'c'){    // catrgorical
                    int cr = rand() % d->domain[c]; // select attribute on which to split
                    cps.push_back(cr);
                    std::vector<int> i0;
                    std::vector<int> i1;
                    for(int i = 0; i < current->indices.size(); i++){
                        int pos = current->indices[i];
                        if(d->data[pos][c] == cr)
                            i1.push_back(pos);
                        else
                            i0.push_back(pos);
                    }
                    indices.push_back(i0);
                    indices.push_back(i1);
                    return T->splitNodeGenerate(indices, cps, toSplit, splitOnID);
                }else{  // integer
                    int procedure = 0;
                    if(tA != 'c' && tA != 'b' && domainShift == 0){
                        int decider = rand() % 3; // 50% prob for split 25% for lin and 25% for quad regression
                        if(decider == 2)
                            procedure = 1;
                        else if(decider == 3)
                            procedure = 2;
                    }
                    if(procedure == 0){         // split
                        int randpos = rand() % ((int) current->indices.size());
                        double cr = d->data[current->indices[randpos]][c];
                        std::vector<int> i0;
                        std::vector<int> i1;
                        cps.push_back(cr);
                        for(int i = 0; i < current->indices.size(); i++){
                            int pos = current->indices[i];
                            if(d->data[pos][c] < cr)
                                i0.push_back(pos);
                            else
                                i1.push_back(pos);
                        }
                        indices.push_back(i0);
                        indices.push_back(i1);
                        return T->splitNodeGenerate(indices, cps, toSplit, splitOnID);
                    }else{   // linear regression
                        int alpha = rand() % domain;
                        int beta = round(((double)(rand() % sqrdom)) / 10.0) + 1;
                        cps.push_back(-1);
                        cps.push_back(alpha);
                        cps.push_back(beta);
                        std::vector<int> i0;
                        for(int i = 0; i < current->indices.size(); i++){
                            int pos = current->indices[i];
                            double v = d->data[pos][c];
                            d->data[pos][a] += alpha;
                            d->data[pos][a] += (beta * v);
                            i0.push_back(pos);
                        }
                        indices.push_back(i0);
                        return T->splitNodeGenerate(indices, cps, toSplit, splitOnID);
                    }/*else{                      // quadratic regression
                        int alpha = rand() % domain;
                        int beta = round(((double)(rand() % sqrdom)) / 10.0);
                        int gamma = round(((double)(rand() % sqrdom)) / 10.0) + 1;
                        cps.push_back(-1);
                        cps.push_back(alpha);
                        cps.push_back(beta);
                        cps.push_back(gamma);
                        std::vector<int> i0;
                        for(int i = 0; i < current->indices.size(); i++){
                            int pos = current->indices[i];
                            double v = d->data[pos][c];
                            d->data[pos][a] += alpha;
                            d->data[pos][a] += (beta * v);
                            d->data[pos][a] += (gamma * v * v);
                            i0.push_back(pos);
                        }
                        indices.push_back(i0);
                        return T->splitNodeGenerate(indices, cps, toSplit, splitOnID);
                    }*/
                }
            }
            currChild++;
        }
    }
    return 0;
}

void generateDataFromLeaves(pns::node* t, Data* d, int a, int domainShift){
    char type = d->getAttributeType(a);
    int categories = 0;
    int domain = 10;
    double prevMean = -1000.0;
    int prevProb = -1;
    if(type == 'c'){
        categories = rand() % 2 + 3;
        d->domain[a] = categories;
    }
    // iterate over leaves
    std::vector<pns::node*> stack;
    stack.push_back(t);
    while(stack.size() > 0){
        pns::node *current = stack.back();
        stack.pop_back();
        for(int i = 0; i < current->children.size(); i++){
            stack.push_back(current->children[i]);
        }
        if(current->isLeaf()){
            if(type == 'b'){
                int rr = 0;
                if(prevProb != -1){
                    rr = (prevProb + 20 + (rand() % 60)) % 100;
                }else{
                    rr = rand() % 60 + 20;
                }
                prevProb = rr;
                for(int i = 0; i < current->indices.size(); i++){
                    int pos = current->indices[i];
                    int newr = rand() % 100;
                    if(newr < rr){
                        d->data[pos][a] = 0.0;
                    }else{
                        d->data[pos][a] = 1.0;
                    }
                }
            }else if(type == 'c'){
                int sum = 0;
                std::vector<int> probs;
                for(int i = 0; i < categories; i++){
                    int rr = rand() % 100;
                    probs.push_back(rr);
                    sum += rr;
                }
                for(int i = 0; i < current->indices.size(); i++){
                    int pos = current->indices[i];
                    int newr = rand() % sum;
                    int probsum = 0;
                    for(int cat = 0; cat < categories; cat++){
                        probsum += probs[cat];
                        if(newr < probsum){
                            d->data[pos][a] = (double)cat;
                            break;
                        }
                    }
                }
            }else{
                double mean = 0.0;
                if(prevMean != -1000.0){
                    int sign = rand() % 2;
                    double shift = ((double)(rand() % domain + 1)) * 0.1;
                    mean = sign == 0 ? prevMean - shift : prevMean + shift;
                }else{
                    mean  = (double) domainShift + rand() % domain;
                }
                prevMean = mean;
                double sd = (double)domainShift + ((double)(rand() % (10 * domain))) / 10.0;
                for(int i = 0; i < current->indices.size(); i++){
                    int pos = current->indices[i];
                    double value = rand_normal(mean, sd);
                    d->data[pos][a] += value;
                }
            }
        }
    }
}
