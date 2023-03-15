#include <cmath>
#include <assert.h>
#include <utility>
#include <set>
#include <map>
#include <algorithm>
#include <limits>

#include "score.h"
////#include "NML_histogram.h"

#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

double binaryLeafModelApproximation(int M, Data *d){
    int p = 10;
    if(M < 1)
    return 0.0;
    if(d->cache[M - 1] != -1)
    return d->cache[M - 1];
    double sum = 1.0;
    double b = 1.0;
    int bound = (int) ceil(2.0 + sqrt(2.0 * M * p * log(10)));
    for(int i = 1; i <= bound; i++){
        b = (M - i + 1) * (b / M);
        sum += b;
    }
    d->cache[M - 1] = sum;
    return d->cache[M - 1];   // changed log to log2 --> also satisfies example from paper
}

double categorialLeafModelApproximation(int M, int K, Data *d){
    double costs = categorialLeafModelApproximationPrecal(M, K, d);
    if(costs <= 0.0)
    return 0.0;
    else
    return log2(costs);
}

double categorialLeafModelApproximationPrecal(int M, int K, Data *d){
    if(K < 1){
        return 0.0;
    }else if(K == 1){
        return 1.0;
    }else{
        double sum = binaryLeafModelApproximation(M, d);
        double old_sum = 1.0;
        if(K > 2){
            for(int j = 3; j <= K; j++){
                double new_sum = sum + (M * old_sum) / ((double)j - 2.0);
                old_sum = sum;
                sum = new_sum;
            }
        }
        return sum;
    }
}


double uniformError(std::vector<double> vals, double resolution){
    double N = (double) vals.size();
    double max = -std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    for(int i = 0; i < vals.size(); i++){
        double v = vals[i];
        if(v < min)
        min = v;
        if(v > max)
        max = v;
    }
    double dataPoint = log2(((max - min) / resolution) + 1);
    return (N * dataPoint);
}

double categorialLeafData(Data* d, pns::node* n){
    if(!n->hasLeafData){
        double score = 0.0;
        if(n->indices.size() == 0)
            return 0.0;
        pns::ad* hm = d->getRestrictedAttributeDistribution(n->attribute, n->indices);
        int N = n->indices.size();
        for(pns::ad::iterator it = hm->begin(); it != hm->end(); it++) {
            if(it->second > 0){
                double temp = (double) it->second / (double) N;
                if(!pns::isZero(temp))
                    score += (-temp * log2(temp));
            }
        }
        delete hm;
        n->leafData = N * score;
        n->hasLeafData = true;
    }
    return n->leafData;
}

double categorialLeafData(Data* d, pns::node* n, pns::ad* hm){
    if(!n->hasLeafData){
        double score = 0.0;
        if(n->indices.size() == 0)
            return 0.0;
        int N = n->indices.size();
        for(pns::ad::iterator it = hm->begin(); it != hm->end(); it++) {
            if(it->second > 0){
                double temp = (double) it->second / (double) N;
                if(!pns::isZero(temp))
                    score += (-temp * log2(temp));
            }
        }
        n->leafData = N * score;
        n->hasLeafData = true;
    }
    return n->leafData;
}

double rvLeafData(Data* d, pns::node* n){
    if(!n->hasLeafData && n->path.size() == 0){
        double N = (double) n->indices.size();
        double score = N * log2(d->domain[n->attribute]);
        n->leafData = score;
        n->hasLeafData = true;
    }else if(!n->hasLeafData){
        double score = 0.0;
        if(n->indices.size() == 0)
            return 0.0;
        double t = 0.0;
        double tt = 0.0;
        double N = (double) n->indices.size();
        std::vector<double> residuals;
        for(int i = 0; i < n->indices.size(); i++){
            int pos = n->indices[i];
            double v = d->getTransformedValue(n, pos);
            residuals.push_back(v);
            t += v;
            tt += (v*v);
        }
        double sse = tt - (((double)1/N) * (t*t));
        if(sse <= 0.0){
            n->leafData = 0.0;
            n->hasLeafData = true;
        }else{
            double sigma = sqrt(sse/N);
            double score = (N/2.0) * log2(sigma*sigma*2.0*M_PI) + ((N/2.0)*(1.0 / log(2.0))) + N * log2(1.0 / d->resolution[n->attribute]);
            score = score < 0.0 ? 0.0 : score;
            double score3 = uniformError(residuals, d->resolution[n->attribute]);
            if(score3 < score){
                score = score3;
            }
            score = score + 1; //uniform or gaussian
            n->leafData = score;
            n->hasLeafData = true;
        }
    }
    return n->leafData;
}

double rvLeafData(Data* d, pns::node* n, double sse){
    if(!n->hasLeafData){
        double score = 0.0;
        if(n->indices.size() == 0)
            return 0.0;
        double N = (double) n->indices.size();
        if(sse <= 0.0){
            n->leafData = 0.0;
            n->hasLeafData = true;
        }else{
            double sigma = sqrt(sse/N);
            double score = (N/2) * log2(sigma*sigma*2.0*M_PI) + ((N/2.0)*(1.0 / log(2.0))) + N * log2(1.0 / d->resolution[n->attribute]);
            score = score < 0.0 ? 0.0 : score;
            n->leafData = score;
            n->hasLeafData = true;
        }
    }
    return n->leafData;
}

double creLeafData(Data* d, pns::node* n){
    if(!n->hasLeafData){
        double score = 0.0;
        if(n->indices.size() == 0)
            return 0.0;
        double N = (double) n->indices.size();
        std::vector<double> values;
        for(int i = 0; i < n->indices.size(); i++){
            int pos = n->indices[i];
            values.push_back(d->getTransformedValue(n, pos));
        }
        std::sort(values.begin(), values.end());
        for(int i = 1; i < values.size(); i++){
            double curr = (values[i] - values[i-1]) * ((double)i/N) * log2((double)i/N);
            score -= curr;
        }
        n->leafData = score;
        n->hasLeafData = true;
    }
    return n->leafData;
}

double CRE(Data* d, std::vector<int>& indices, int a){
    double score = 0.0;
    if(indices.size() == 0)
        return 0.0;
    double N = (double) indices.size();
    std::vector<double> values;
    for(int i = 0; i < indices.size(); i++){
        int pos = indices[i];
        double v = d->data[pos][a];
        values.push_back(d->data[pos][a]);
    }
    std::sort(values.begin(), values.end());
    for(int i = 1; i < values.size(); i++){
        double curr = (values[i] - values[i-1]) * ((double)i/N) * log2((double)i/N);
        score -= curr;
    }
    double lgscore = log2((score + 1) / d->resolution[a]);
    return lgscore;
}

double rvLeafModel(Data* d, pns::node* n){
    if(n->path.size() == 0){
        return 0.0;
    }else{
        return 2.0 * log2(d->domain[n->attribute]);
    }
}

double treeScore(Data* d, pns::node* t, int a){
    double score = 0.0;
    std::vector<pns::node*> stack;
    stack.push_back(t);
    double constantInternal = 2.0 + log2(d->numA);
    while(stack.size() > 0){
        pns::node *current = stack.back();
        stack.pop_back();
        for(int i = 0; i < current->children.size(); i++){
            stack.push_back(current->children[i]);
        }
        char type = d->getAttributeType(current->attribute);
        if(current->isLeaf()){
            score += leafCosts(d, current);
        }else{
            if(type == 'c'){
                if(current->multiwaySplit){
                    score += constantInternal;
                }else{
                    score += constantInternal + log2(d->domain[current->attribute]);
                }
            }else if(current->isRegressor){
                score += constantInternal + regressionModelCosts(d, current);
            }else if(type == 'b'){
                score += constantInternal;
            }else{
                int n0 = pns::min((int)current->children[0]->indices.size(), (int)current->children[1]->indices.size());
                score += constantInternal + log2(d->domain[current->attribute] - 1);///*CRE(d, current->indices, a);/*/getSplitPointCosts(d->rows/*/current->indices.size()*/, n0); //+ log2(d->domain[current->attribute] - 1);
            }
        }
    }
    return score;
}

double forestScore(Data* d, std::vector<Tree*> forest){
    double score = 0.0;
    for(int i = 0; i < forest.size(); i++){
        Tree* tt = forest[i];
        double sc = treeScore(d, tt->getTree(), tt->getPrimeAttribute());
        score += sc;
    }
    return score;
}

double leafData(Data* d, pns::node* n){
    char t = d->getAttributeType(n->attribute);
    if(t == 'b' || t == 'c'){
        return categorialLeafData(d,n);
    }else{
        return rvLeafData(d,n);
    }
}

double leafModelApproximation(Data *d, pns::node* n){
    char t = d->getAttributeType(n->attribute);
    if(t == 'b' || t == 'c')
        return categorialLeafModelApproximation(n->indices.size(), d->domain[n->attribute], d);
    else
        return rvLeafModel(d,n);
}

double leafCosts(Data *d, pns::node* n){
    return 1.0 + leafModelApproximation(d,n) + leafData(d,n);
}

pns::IndCut* splitRegression(Data* d, pns::node* n, int cand){
    // sort indices according to candidate
    std::vector<pns::IndValPair> mapping;
    for(int i = 0; i < n->indices.size(); i++){
        pns::IndValPair p;
        p.pos = n->indices[i];
        p.value = d->data[p.pos][cand];
        mapping.push_back(p);
    }
    std::sort(mapping.begin(), mapping.end(), pns::by_value());
    // find efficiently best position and costs
    pns::IndValPair best;
    best = regressionScoring(d, n, mapping);
    // calculate score
    std::vector<int> zeroD;
    std::vector<int> oneD;
    double best_cp = 0.0;
    for(int i = 0; i < mapping.size(); i++){
        if(i <= best.pos)
            oneD.push_back(mapping[i].pos);
        else{
            if(i == (best.pos + 1))
                best_cp = mapping[i].value;
            zeroD.push_back(mapping[i].pos);
        }
    }
    double newCosts = splitLeafCosts(zeroD, oneD, d, n);
    best.value = newCosts;
    double oldCosts = leafCosts(d, n);
    double finalCosts = newCosts + 2.0 - oldCosts + log2(d->numA) + log2(d->domain[cand] - 1);
    
    pns::IndValPair pairMultiway = regressionScoringMultiway(d, n, mapping);
    double scoreMultiway = pairMultiway.value + 2.0 - oldCosts + log2(d->numA) + log2(10);
    
    bool applyRegression = false;
    pns::Regressor* regressor = NULL;
    char type = d->getAttributeType(n->attribute);
    // check if regressor can be trained
    if(type != 'b' && type != 'c' && d->domain[cand] > 2){
        // calculate regressor and sse for quadratic regression
        pns::Regressor* regressorQ = calculateQuadraticRegressor(d, n, mapping);
        if(regressorQ){
            pns::node* curr = copyNode(n, n->indices);
            // get leaf cost (gaussian point model wrt. applied regression)
            double regressionLeaf = 1.0 + rvLeafModel(d, n);
            regressionLeaf += rvLeafData(d, curr, regressorQ->sse);
            delete curr;
            // costs for regression internal node
            double reginternal = 2.0 + log2(d->numA) + regressionModelCosts(d, regressorQ);
            // delta
            double regCosts = regressionLeaf + reginternal - oldCosts;
            if(regCosts < finalCosts){
                finalCosts = regCosts;
                applyRegression = true;
                regressor = regressorQ;
            }
        }
        // calculate regressor and sse for linear regression
        pns::Regressor* regressorL = calculateRegressor(d, n, mapping);
        if(regressorL){
            pns::node* curr = copyNode(n, n->indices);
            // get leaf cost (gaussian point model wrt. applied regression)
            double regressionLeaf = 1.0 + rvLeafModel(d, n);
            regressionLeaf += rvLeafData(d, curr, regressorL->sse);
            delete curr;
            // costs for regression internal node
            double reginternal = 2.0 + log2(d->numA) + regressionModelCosts(d, regressorL);
            // delta
            double regCosts = regressionLeaf + reginternal - oldCosts;
            if(regCosts < finalCosts){
                finalCosts = regCosts;
                applyRegression = true;
                regressor = regressorL;
            }
        }
    }
    std::vector<double> cps;
    std::vector<std::vector<int> > ret;
    if(scoreMultiway < finalCosts){ // multiway split
        finalCosts = scoreMultiway;
        double K = pairMultiway.pos;
        cps.push_back(K);
        std::map<double, int> counts;
        for(int i = 0; i < mapping.size(); i++){
            counts[mapping[i].value]++;
        }
        // initialize bins
        int index = 1;
        std::vector<int> small_ones;
        ret.push_back(small_ones);
        std::map<double, int> pointer;
        std::map<double, int>::iterator it;
        for(it = counts.begin(); it != counts.end(); it++){
            double key = it->first;
            double val = it->second;
            if(val < K){
                pointer[key] = 0;
            }else{
                pointer[key] = index;
                std::vector<int> dummy;
                ret.push_back(dummy);
                index++;
            }
        }
        // fill bins
        for(int i = 0; i < mapping.size(); i++){
            ret[pointer[mapping[i].value]].push_back(mapping[i].pos);
        }
    }else if(applyRegression){
        std::vector<int> one;
        for(int i = 0; i < mapping.size(); i++){
            one.push_back(mapping[i].pos);
        }
        cps.push_back(-1.0);    //indicator for regressor
        cps.push_back(regressor->alpha);
        cps.push_back(regressor->beta);
        if(regressor->hasGamma)
            cps.push_back(regressor->gamma);
        delete regressor;
        ret.push_back(one);
    }else{
        std::vector<int> zero;
        std::vector<int> one;
        double best_cp = 0.0;
        for(int i = 0; i < mapping.size(); i++){
            if(i <= best.pos)
                one.push_back(mapping[i].pos);
            else{
                if(i == (best.pos + 1))
                    best_cp = mapping[i].value;
                zero.push_back(mapping[i].pos);
            }
        }
        cps.push_back(best_cp); // le
        ret.push_back(one);
        ret.push_back(zero);
    }
    pns::IndCut* dd = new pns::IndCut();
    dd->first = ret;
    dd->second.first = cps;
    dd->second.second = finalCosts;
    return dd;
}

double splitCategorialComplete(Data* d,pns::node* n, std::vector<pns::IndValPair>& mapping){
    pns::IndValPair pair;
    double currentValue = 0.0;
    bool valueSet = false;
    double costs = 0.0;
    std::vector<int> current;
    int count = 0;
    for(int i = 0; i < mapping.size(); i++){
        if(!valueSet){
            currentValue = mapping[i].value;
            current.push_back(mapping[i].pos);
            valueSet = true;
        }else{
            if(mapping[i].value != currentValue){
                if(current.size() > 1)
                    count++;
                pns::node* n0 = copyNode(n, current);
                double currentCosts = leafCosts(d,n0);
                costs += currentCosts;
                std::vector<int> dummy;
                currentValue = mapping[i].value;
                current = dummy;
                current.push_back(mapping[i].pos);
            }else{
                current.push_back(mapping[i].pos);
            }
        }
        if((i+1) == mapping.size() && current.size() > 1){
            pns::node* n0 = copyNode(n, current);
            double currentCosts = leafCosts(d,n0);
            costs += currentCosts;
            count++;
        }
    }
    if(count > 2){
        return costs;
    }else{
        return std::numeric_limits<double>::max();;
    }
}

pns::IndCut* splitCategorial(Data* d, pns::node* n, int cand){
    std::vector<pns::IndValPair> mapping;
    for(int i = 0; i < n->indices.size(); i++){
        pns::IndValPair p;
        p.pos = n->indices[i];
        p.value = d->data[p.pos][cand];
        mapping.push_back(p);
    }
    std::sort(mapping.begin(), mapping.end(), pns::by_value());
    double multiwaySplitScore = splitCategorialComplete(d,n,mapping);
    multiwaySplitScore += + 2.0 + log2(d->numA);
    char typeAtt = d->getAttributeType(n->attribute);
    std::vector<int> zero;
    std::vector<int> one;
    std::vector<int> old = n->indices;
    assert(old.size() > 0);
    double cat1 = d->data[old[0]][cand];
    std::set<double> categs;
    categs.insert(cat1);
    // iter first time to get categories at the same time
    for(int i = 0; i < old.size(); i++){
        double curr = d->data[old[i]][cand];
        if(curr == cat1){
            one.push_back(old[i]);
        }else{
            zero.push_back(old[i]);
        }
        categs.insert(curr);
    }
    double best = cat1;
    double bestScore = splitLeafCosts(zero, one, d, n);
    typedef std::set<double>::iterator cit;
    cit pos = categs.find(cat1);
    categs.erase(pos);
    while(categs.size() > 0){
        cit cur = categs.begin();
        std::vector<int> zeroD;
        std::vector<int> oneD;
        double cat = *cur;
        for(int i = 0; i < old.size(); i++){
            if(d->data[old[i]][cand] == cat){
                oneD.push_back(old[i]);
            }else{
                zeroD.push_back(old[i]);
            }
        }
        double score = splitLeafCosts(zeroD, oneD, d, n);
        if(score < bestScore){
            bestScore = score;
            zero = zeroD;
            one = oneD;
            best = cat;
        }
        categs.erase(cur);
    }
    double newCosts = bestScore + 2.0 + log2(d->numA) + log2(d->domain[cand]);
    double oldCosts = leafCosts(d, n);
    pns::IndCut* dd = new pns::IndCut();
    std::vector<double> cps;
    std::vector<std::vector<int> > ret;
    if(multiwaySplitScore < newCosts){
        newCosts = multiwaySplitScore;
        // add all indices
        double currentValue = 0.0;
        bool valueSet = false;
        std::vector<int> current;
        for(int i = 0; i < mapping.size(); i++){
            if(!valueSet){
                currentValue = mapping[i].value;
                current.push_back(mapping[i].pos);
                valueSet = true;
            }else{
                if(mapping[i].value != currentValue){
                    ret.push_back(current);
                    std::vector<int> dummy;
                    currentValue = mapping[i].value;
                    current = dummy;
                    current.push_back(mapping[i].pos);
                }else{
                    current.push_back(mapping[i].pos);
                }
            }
            if((i+1) == mapping.size() && current.size() > 1){
                ret.push_back(current);
            }
        }
        cps.push_back(-1.0);
    }else{
        ret.push_back(zero);
        ret.push_back(one);
        cps.push_back(best);
    }
    dd->first = ret;
    dd->second.first = cps;
    // calculate score
    dd->second.second = newCosts - oldCosts;
    return dd;
}

pns::IndCut* splitBinary(Data* d, pns::node* n, int cand){
    std::vector<int> zero;
    std::vector<int> one;
    std::vector<int> old = n->indices;
    for(int i = 0; i < old.size(); i++){
        if(d->data[old[i]][cand] == 0){
            zero.push_back(old[i]);
        }else{
            one.push_back(old[i]);
        }
    }
    std::vector<std::vector<int> > ret;
    ret.push_back(zero);
    ret.push_back(one);
    std::vector<double> cps;
    cps.push_back(1);
    pns::IndCut* dd = new pns::IndCut();
    dd->first = ret;
    dd->second.first = cps;
    double newCosts = splitLeafCosts(zero, one, d, n);
    double oldCosts = leafCosts(d, n);
    dd->second.second = newCosts + 2.0 + log2(d->numA) - oldCosts;
    return dd;
}

double splitLeafData(std::vector<int>& zero, std::vector<int>& one, Data* d,pns::node* n){
    pns::node* n0 = copyNode(n, one);
    pns::node* n1 = copyNode(n, zero);
    double result = leafData(d, n0) + leafData(d, n1);
    delete n0;
    delete n1;
    return result;
}

double splitLeafModel(std::vector<int>& zero, std::vector<int>& one, Data* d, pns::node* n){
    pns::node* n0 = copyNode(n, one);
    pns::node* n1 = copyNode(n, zero);
    double result = leafModelApproximation(d, n0) + leafModelApproximation(d, n1);
    delete n0;
    delete n1;
    return result;
}

double splitLeafCosts(std::vector<int>& zero, std::vector<int>& one, Data* d, pns::node* n){
    pns::node* n0 = copyNode(n, one);
    pns::node* n1 = copyNode(n, zero);
    double result = leafCosts(d, n0) + leafCosts(d, n1);
    delete n0;
    delete n1;
    return result;
}

double splitLeafCosts(std::vector<int>& zero,std::vector<int>& one,Data* d,pns::ad* ad0,pns::ad* ad1, pns::node* n){
    pns::node* n0 = copyNode(n, one);
    pns::node* n1 = copyNode(n, zero);
    double result = 2 + leafModelApproximation(d, n0) + leafModelApproximation(d, n1) + categorialLeafData(d,n0,ad1) + categorialLeafData(d,n1,ad0);
    delete n0;
    delete n1;
    return result;
}

double splitLeafCosts(std::vector<int>& zero,std::vector<int>& one,Data* d, double sse0,double sse1, pns::node* n){
    pns::node* n0 = copyNode(n, one);
    pns::node* n1 = copyNode(n, zero);
    double result = 2 + leafModelApproximation(d, n0) + leafModelApproximation(d, n1);
    result += rvLeafData(d,n0,sse1) + rvLeafData(d,n1,sse0);
    delete n0;
    delete n1;
    return result;
}

pns::IndCut* splitNode(Data* d, pns::node* n , int cand){
    char t = d->getAttributeType(cand);
    if(t == 'b')
        return splitBinary(d,n,cand);
    if(t == 'c')
        return splitCategorial(d,n,cand);
    else{
        return splitRegression(d,n,cand);
    }
}


pns::IndValPair regressionScoring(Data* d,pns::node* a,std::vector<pns::IndValPair>& mapping){
    char t = d->getAttributeType(a->attribute);
    if(t == 'b' || t == 'c')
        return regressionScoringCB(d,a,mapping);
    else
        return regressionScoringReg(d,a,mapping);
}

pns::IndValPair regressionScoringReg(Data* d,pns::node* n,std::vector<pns::IndValPair>& mapping){
    int a = n->attribute;
    std::vector<int> zero;
    std::vector<int> one;
    double best_cp = 0.0;
    for(int i = 0; i < mapping.size(); i++){
        zero.push_back(mapping[i].pos);
    }
    // precalculate scores for fast sse calculation
    std::vector<double> cy;
    std::vector<double> cyy;
    double t = 0.0;
    double tt = 0.0;
    cy.push_back(t);
    cyy.push_back(tt);
    double N = (double) zero.size();
    for(int i = 0; i < zero.size(); i++){
        int pos = zero[i];
        double v = d->getTransformedValue(n, pos);
        double myx = mapping[i].value;
        t += v;
        tt += (v * v);
        cy.push_back(t);
        cyy.push_back(tt);
    }
    pns::IndValPair best;
    best.pos = 0; // position up to which elements are included to one
    best.end = -1;
    best.value = std::numeric_limits<double>::max(); // only negative values are of benefit
    int iN = mapping.size();
    for(int k = 0; k < iN - 1; k++){
        one.push_back(mapping[k].pos);
        zero.erase(zero.begin());
        if((k+2) == mapping.size() || mapping[k].value != mapping[k+1].value){
            // calculate sse from 0 to i (everything +1 cause cx and cxx have one more element)
            int i = 1;
            int j = k+1;
            double yyij = cyy[j] - cyy[i-1];
            double yij = cy[j] - cy[i-1];
            double elems = (double)j-i+1;
            double sse1 = (yyij) - (1.0/(elems)) * (yij * yij);
            // calculate sse from i+1 to end
            i = k+2;
            j = mapping.size();
            yyij = cyy[j] - cyy[i-1];
            yij = cy[j] - cy[i-1];
            elems = (double)j-i+1;
            double sse0 = (yyij) - (1.0/(elems)) * (yij * yij);
            double currentCosts = splitLeafCosts(zero, one, d, sse0, sse1, n); //+ getSplitPointCosts(d->rows/*/iN*/, pns::min(k+1, mapping.size() - (k+1)));
            if(currentCosts < best.value){
                best.value = currentCosts;
                best.pos = k;
            }
        }
    }
    return best;
}

pns::IndValPair regressionScoringMultiway(Data* d,pns::node* n,std::vector<pns::IndValPair>& mapping){
    int a = n->attribute;
    pns::IndValPair best;
    best.pos = 0; // position up to which elements are included to one
    best.end = -1;
    best.value = std::numeric_limits<double>::max(); // only negative values are of benefit
    // count duplicates
    std::map<double, int> counts;
    for(int i = 0; i < mapping.size(); i++){
        counts[mapping[i].value]++;
    }
    //iterate over possible counts
    for(int k = 2; k <= 10; k++){
        // initialize bins
        int index = 1;
        std::vector<std::vector<int>> positions;
        std::vector<int> small_ones;
        positions.push_back(small_ones);
        int need_small_ones = 0;
        std::map<double, int> pointer;
        std::map<double, int>::iterator it;
        for(it = counts.begin(); it != counts.end(); it++){
            double key = it->first;
            double val = it->second;
            if(val < k){
                pointer[key] = 0;
                need_small_ones++;
            }else{
                pointer[key] = index;
                std::vector<int> dummy;
                positions.push_back(dummy);
                index++;
            }
        }
        int contIndices = need_small_ones >= 2 ? positions.size() : positions.size() - 1;
        if(contIndices < 2){
            break;
        }
        // fill bins
        for(int i = 0; i < mapping.size(); i++){
            if(need_small_ones < 2){
                if(pointer[mapping[i].value] == 0){
                    positions[1].push_back(mapping[i].pos);
                }else{
                    positions[pointer[mapping[i].value]].push_back(mapping[i].pos);
                }
            }else{
                positions[pointer[mapping[i].value]].push_back(mapping[i].pos);
            }
        }
        // calculate costs
        double costs = 0;
        for(int i = 0; i < positions.size(); i++){
            if(positions[i].size() == 1){
            }
            if(positions[i].size() < 2)
                continue;
            pns::node* n1 = copyNode(n, positions[i]);
            if(i != 0)
                n1->duplicates = true;
            costs += leafCosts(d, n1);
        }
        if(costs < best.value){
            best.value = costs;
            best.pos = k;
        }
    }
    return best;
}

pns::IndValPair regressionScoringInterval(Data* d,pns::node* n,std::vector<pns::IndValPair>& mapping){
    int a = n->attribute;
    std::vector<int> zero;
    std::vector<int> one;
    std::vector<double> vals;
    double best_cp = 0.0;
    for(int i = 0; i < mapping.size(); i++){
        int pos = mapping[i].pos;
        zero.push_back(pos);
        double v = d->getTransformedValue(n, pos);
        vals.push_back(v);
        
    }
    pns::IndValPair best;
    best.pos = 0; // position up to which elements are included to one
    best.end = -1;
    best.value = std::numeric_limits<double>::max(); // only negative values are of benefit
    int iN = mapping.size();
    for(int i = 1; i < iN - 3; i++){
        if(vals[i] == vals[i-1])
            continue;
        std::vector<int> one_buff;
        std::vector<int> zero_buff;
        for(int j = i+2; j < iN; j++){ // index j is excluded
            if(vals[j] == vals[j-1])
                continue;
            // calculate current sses
            double zt = 0.0;
            double ztt = 0.0;
            double ot = 0.0;
            double ott = 0.0;
            for(int k = 0; k < iN; k++){
                double v = vals[k];
                double vv = v * v;
                if(k >= i && k < j){
                    one_buff.push_back(zero[k]);
                    ot += v;
                    ott += vv;
                }else{
                    zero_buff.push_back(zero[k]);
                    zt += v;
                    ztt += vv;
                }
            }
            double sse0 = ztt - ((1.0/(double)zero_buff.size()) * (zt*zt));
            double sse1 = ott - ((1.0/(double)one_buff.size()) * (ot*ot));
            double currentCosts = splitLeafCosts(zero_buff, one_buff, d, sse0, sse1, n); //+ getSplitPointCosts(d->rows/*/iN*/, pns::min(k+1, mapping.size() - (k+1)));
            if(currentCosts < best.value){
                best.value = currentCosts;
                best.pos = i;
                best.end = j;
            }
        }
    }
    return best;
}

pns::IndValPair regressionScoringCB(Data* d,pns::node* n,std::vector<pns::IndValPair>& mapping){
    int a = n->attribute;
    std::vector<int> zero;
    std::vector<int> one;
    double best_cp = 0.0;
    for(int i = 0; i < mapping.size(); i++){
        zero.push_back(mapping[i].pos);
    }
    // track counts for current position
    pns::ad* ad0 = d->getRestrictedAttributeDistribution(a, zero);
    pns::ad* ad1 = new pns::ad();
    
    pns::IndValPair best;
    best.pos = 0; // position up to which elements are included to one
    best.end = -1;
    best.value = std::numeric_limits<double>::max(); // only negative values are of benefit
    int N = mapping.size();
    for(int i = 0; i < N - 1; i++){
        // change counts
        double cval = d->data[mapping[i].pos][a];
        (*ad0)[cval] -= 1;
        (*ad1)[cval] += 1;
        one.push_back(mapping[i].pos);
        zero.erase(zero.begin());
        if((i+2) == mapping.size() || mapping[i].value != mapping[i+1].value){
            double currentCosts = splitLeafCosts(zero, one, d, ad0, ad1, n); //+ getSplitPointCosts(d->rows/*/N*/, pns::min(i+1, mapping.size() - (i+1)));
            if(currentCosts < best.value){
                best.value = currentCosts;
                best.pos = i;
            }
        }
    }
    delete ad0;
    delete ad1;
    return best;
}

pns::IndValPair regressionScoringRegNML(Data* d,pns::node* n,std::vector<pns::IndValPair>& mapping){
    std::vector<int> zero;
    std::vector<int> one;
    int a = n->attribute;
    double best_cp = 0.0;
    for(int i = 0; i < mapping.size(); i++){
        zero.push_back(mapping[i].pos);
    }
    // fix step size
    int stepSize = (int) ceil(mapping.size() / pow((double) mapping.size(), d->binBound));
    stepSize = stepSize >= 1 ? stepSize : 1;
    pns::IndValPair best;
    best.pos = 0; // position up to which elements are included to one
    best.end = -1;
    best.value = std::numeric_limits<double>::max(); // only negative values are of benefit
    for(int i = stepSize; i < mapping.size() - 1; i += stepSize){
        for(int k = i - stepSize; k < i; k++){
            one.push_back(mapping[k].pos);
            zero.erase(zero.begin());
        }
        double currentCosts = splitLeafCosts(zero, one, d, n); //+ getSplitPointCosts(d->rows/*/iN*/, pns::min(i, mapping.size() - (i)));
        if(currentCosts < best.value){
            best.value = currentCosts;
            best.pos = (i - 1);
        }
    }
    return best;
}

pns::Regressor* calculateRegressor(Data* d,pns::node* n,std::vector<pns::IndValPair>& mapping){
    int a = n->attribute;
    double ty = 0.0;
    double xy = 0.0;
    double tx = 0.0;
    double xx = 0.0;
    double N = (double) mapping.size();
    for(int i = 0; i < mapping.size(); i++){
        int pos = mapping[i].pos;
        double v = d->getTransformedValue(n, pos);
        double myx = mapping[i].value;
        ty += v;
        tx += myx;
        xx += myx * myx;
        xy += myx * v;
    }
    double tyN = ty / N;
    double xyN = xy / N;
    double xxN = xx / N;
    double txN = tx / N;
    double nom = (xyN - (txN * tyN));
    double denom = (xxN - (txN * txN));
    if(pns::isZero(nom) || pns::isZero(denom))
        return NULL;
    double beta = nom / denom;
    double alpha = tyN - ((tx * beta)/N);
    double t = 0.0;
    double tt = 0.0;
    for(int i = 0; i < mapping.size(); i++){
        int pos = mapping[i].pos;
        double v = d->getTransformedValue(n, pos);
        double yhead = alpha + (mapping[i].value * beta);
        double err = v - yhead;
        t += err;
        tt += (err * err);
    }
    double sse = tt - ((1.0/N) * (t*t));
    pns::Regressor* reg = new pns::Regressor();
    reg->beta = beta;
    reg->alpha = alpha;
    reg->hasGamma = false;
    reg->sse = sse;
    return reg;
}

pns::Regressor* calculateQuadraticRegressor(Data* d,pns::node* n,std::vector<pns::IndValPair>& mapping){
    int a = n->attribute;
    double ty = 0.0;
    double xy = 0.0;
    double tx = 0.0;
    double xx = 0.0;
    double x3 = 0.0;
    double x4 = 0.0;
    double x2y = 0.0;
    double N = (double) mapping.size();
    for(int i = 0; i < mapping.size(); i++){
        int pos = mapping[i].pos;
        double v = d->getTransformedValue(n, pos);
        double myx = mapping[i].value;
        double myx2 = myx * myx;
        ty += v;
        tx += myx;
        xx += myx2;
        x3 += myx * myx2;
        x4 += myx2 * myx2;
        xy += myx * v;
        x2y += myx2 * v;
    }
    double tyN = ty / N;
    double xyN = xy / N;
    double xxN = xx / N;
    double txN = tx / N;
    double x3N = x3 / N;
    double x4N = x4 / N;
    double x2yN = x2y / N;
    double Sxx = xxN - (txN * txN);
    double Sxy = xyN - (txN * tyN);
    double Sxx2 = x3N - (txN * xxN);
    double Sx2x2 = x4N - (xxN * xxN);
    double Sx2y = x2yN - (xxN * tyN);
    double denom = Sxx * Sx2x2 - (Sxx2 * Sxx2);
    if(pns::isZero(denom))
        return NULL;
    double beta = (Sxy * Sx2x2 - Sx2y * Sxx2) / denom;
    double gamma = (Sx2y * Sxx - Sxy * Sxx2) / denom;
    if(pns::isZero(gamma) && pns::isZero(beta))
        return NULL;
    double alpha = tyN - beta * txN - gamma * xxN;
    double t = 0.0;
    double tt = 0.0;
    for(int i = 0; i < mapping.size(); i++){
        int pos = mapping[i].pos;
        double v = d->getTransformedValue(n, pos);
        double mv = mapping[i].value;
        double yhead = alpha + (mv * beta) + ((mv * mv) * gamma);
        double err = v - yhead;
        t += err;
        tt += (err * err);
    }
    double sse = tt - ((1.0/N) * (t*t));
    pns::Regressor* reg = new pns::Regressor();
    reg->beta = beta;
    reg->alpha = alpha;
    reg->hasGamma = true;
    reg->gamma = gamma;
    reg->sse = sse;
    return reg;
}

pns::node* copyNode(pns::node* n, std::vector<int>& indices){
    pns::node* curr = new pns::node();
    curr->attribute = n->attribute;
    curr->indices = indices;
    curr->hasLeafData = false;
    curr->alpha = n->alpha;
    curr->beta = n->beta;
    curr->gamma = n->gamma;
    curr->beta_ref = n->beta_ref;
    curr->path = n->path;
    curr->path.push_back(n->attribute); // dummy
    return curr;
}

pns::node* copyNode(pns::node* n, std::vector<int>& indices, pns::Regressor* reg, int cand){
    pns::node* curr = new pns::node();
    curr->attribute = n->attribute;
    curr->indices = indices;
    curr->hasLeafData = false;
    curr->alpha = n->alpha + reg->alpha;
    curr->beta = n->beta;
    curr->beta.push_back(reg->beta);
    curr->gamma = n->gamma;
    if(reg->hasGamma){
        curr->gamma.push_back(reg->gamma);
    }else{
        curr->gamma.push_back(0.0);
    }
    curr->beta_ref = n->beta_ref;
    curr->beta_ref.push_back(cand);
    curr->path = n->path;
    curr->path.push_back(n->attribute); // dummy
    return curr;
}

double getSplitPointCosts(int in, int in0){
    double n = (double) in;
    double n0 = (double) in0;
    return -log2(n0 / n);
}

double paramCost(double val){
    val = pns::absD(val);
    double prec = 1;
    if(val != 0){
        while(val < 1000){
            val *= 10;
            prec += 1;
            if(prec >= 1000)
                break;
        }
    }
    return 1 + pns::log2N(prec) + pns::log2N(val);
}

double regressionModelCosts(Data* d, pns::node* n){
    double costs = 0.0;
    for(int i = 1; i < n->cut_points.size(); i++){
        costs += paramCost(n->cut_points[i]);
    }
    return costs;
}

double regressionModelCosts(Data* d, pns::Regressor* r){
    double costs = 0.0;
    costs += paramCost(r->alpha);
    costs += paramCost(r->beta);
    if(r->hasGamma){
        costs += paramCost(r->gamma);
    }
    return costs;
}
