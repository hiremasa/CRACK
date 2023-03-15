#ifndef ORIGI_H
#define ORIGI_H

#include <vector>
#include <utility>

#include "forest.h"
#include "score.h"
#include "Wall_Time.h"

namespace ERGO {
    
    /**
     * whichOne is  0   compress L(X)
     *              1   compress L(Y)
     *              2   compress L(X|Y) (Y is right)
     *              3   compress L(Y|X) (X is left)
     **/
    pns::Result inline runGreedyAlgorithm(Data *data, int precision, int whichOne){
        double t1 = getTime();
        // create initial forest
        Forest f(data);
        // create trivial trees
        int start = whichOne % 2 != 0 ? data->xColsSmallerThan : 0;
        int end = whichOne % 2 == 0 ? data->xColsSmallerThan : data->numA;
        fprintf(stdout, "start: %d, end: %d\n", start, end);
        for(int i = start; i < end; i++){
            if(whichOne == 0){  // set up candidates for L(X) + L(Y)
                f.createTrivialTree(i, 0);
            }else if(whichOne == 1){
                f.createTrivialTree(i, 1);
            }else if(whichOne == 2){
                f.createTrivialTree(i, 1);
            }else if(whichOne == 3){
                f.createTrivialTree(i, 0);
            }else if(whichOne == 4){
                f.createTrivialTree(i, 1);
            }else if(whichOne == 5){
                f.createTrivialTree(i, 0);
            }else{
                f.createTrivialTree(i, 2);
            }
        }
        double initialScore = 0.0;
        if(data->nci){
            initialScore = f.size();
        }else{
            initialScore = f.getForestScore();
        }
        double sumChange = 0.0; //
        int numOfSplits = 0;
        std::cout << "------> Initial score: " << std::to_string(initialScore) << "\n";
        bool change = true;
        bool secondDone = false;
        while(change){
            change = false;
            pns::candidate cand;
            cand.init();
            // search for best candidate split among all trees
            for(int i = start; i < end; i++){
                pns::candidate currentCand = f.findCandidateSplit(i);
                if(currentCand.splitOn == -1)
                    continue;
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
            }else{
                if(!secondDone){
                    secondDone = true;
                    if(whichOne == 4){
                        fprintf(stdout, "===============\n");
                        f.addOtherOnes(1);
                        change = true;
                    }else if(whichOne == 5){
                        fprintf(stdout, "===============\n");
                        f.addOtherOnes(0);
                        change = true;
                    }
                }
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
