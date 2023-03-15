#ifndef PACK_H
#define PACK_H

#include <vector>
#include <utility>

#include "forest.h"
#include "score.h"
#include "Wall_Time.h"

namespace RACK {
    
    /**
     * whichOne is  -1  if no causal inference is done
     *              0   compress each by themselves (L(X) + L(Y))
     *              1   compress L(X) + L(Y|X) (X is left)
     *              2   compress L(Y) + L(X|Y) (Y is right)
     *              3   compress attribute x given all other attributes
     **/
    pns::Result inline runGreedyAlgorithm(Data *data, int precision, int whichOne){
        double t1 = getTime();
        // create initial forest
        Forest f(data);
        // create trivial trees
        int start = whichOne == 3 ? (data->xColsSmallerThan - 1) : 0;
        int end = whichOne == 3 ? data->xColsSmallerThan : data->numA;
        for(int i = 0; i < data->numA; i++){
            if(whichOne == 0){  // set up candidates for L(X) + L(Y)
                if(i < data->xColsSmallerThan){
                    f.createTrivialTree(i, 0);
                }else{
                    f.createTrivialTree(i, 1);
                }
            }else if(whichOne == 1){
                if(i < data->xColsSmallerThan){
                    f.createTrivialTree(i, 0);
                }else{
                    f.createTrivialTree(i, 2);
                }
            }else if(whichOne == 2){
                if(i < data->xColsSmallerThan){
                    f.createTrivialTree(i, 2);
                }else{
                    f.createTrivialTree(i, 1);
                }
            }else{
                f.createTrivialTree(i, whichOne);
            }
        }
        double initialScore = f.getForestScore();
        double sumChange = 0.0; //
        int numOfSplits = 0;
        std::cout << "------> Initial score: " << std::to_string(initialScore) << "\n";
        bool change = true;
        while(change){
            change = false;
            pns::candidate cand;
            cand.init();
            // search for best candidate split among all trees
            for(int i = start; i < end; i++){
                pns::candidate currentCand = f.findCandidateSplit(i);
                if(currentCand.splitOn == -1){
                    continue;
                }
                if(currentCand.change <= cand.change){
                    change = true;
                    cand.set(currentCand);
                }
            }
            // perform split
            if(change){
                sumChange += cand.change; //
                if(cand.ic->first.size() > 2){
                    std::cout << "Multiway split at " << cand.toSplit->attribute << " with " << cand.splitOn << "; change: " << cand.change << std::endl;
                }else if(cand.ic->second.first[0] == -1.0 && cand.ic->second.first.size() == 3)
                    std::cout << "Regress " << cand.toSplit->attribute << " with " << cand.splitOn << "; change: " << cand.change << std::endl;
                else if(cand.ic->second.first[0] == -1.0 && cand.ic->second.first.size() == 4)
                    std::cout << "Quadratic regression on " << cand.toSplit->attribute << " with " << cand.splitOn << "; change: " << cand.change << std::endl;
                else
                    std::cout << "Split at " << cand.toSplit->attribute << " with " << cand.splitOn << "; change: " << cand.change << std::endl;
                f.addDAG(cand.splitOn, cand.toSplit->attribute);
                f.split(cand);
                f.removeCandidates();
                numOfSplits++;
            }
        }
        double changeScore = initialScore + sumChange;
        double finalScore = changeScore;
        fprintf(stdout, "Performed %d compression steps.\n", numOfSplits);
        std::string comp = "Compression: " + std::to_string(finalScore) + " / " + std::to_string(initialScore) + " = " + std::to_string(finalScore/initialScore);
        std::cout << "------> " << comp << "\n";
        double execution_time = getTime() - t1;
        std::cout << "Run time: " << execution_time << " seconds (wall time)." << std::endl;
        // write output
        pns::Result r;
        std::vector< pns::dot_file* > files = f.forestToDot(precision);
        r.set(files, initialScore, finalScore, execution_time);
        return r;
    }
}

#endif
