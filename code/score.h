#ifndef SCORE_H
#define SCORE_H

#include <vector>
#include <map>

#include "tree.h" // sould also include defs and data

/**
 * Methods to calculate the model costs of a leaf
 **/
//double binaryLeafModelApproximation(int,Data*);
//double categorialLeafModelApproximation(int, int, Data*);
//double categorialLeafModelApproximationPrecal(int, int, Data*);
// those are moved to the NML file
double rvLeafModel(Data* d, pns::node* n);

/**
 * Uniform error
 **/
double uniformError(std::vector<double>,double);

/**
 * Model costs for regressor
 **/
double regressionModelCosts(Data*, pns::Regressor*);
double regressionModelCosts(Data*, pns::node*);

/**
 * Methods to caluclate the data costs of a leaf
 **/
double categorialLeafData(Data*, pns::node*);
double categorialLeafData(Data*, pns::node*, pns::ad*);
double rvLeafData(Data*, pns::node*);
double rvLeafData(Data*,pns::node*,double);
double creLeafData(Data*, pns::node*);
double CRE(Data*, std::vector<int>&, int);


/**
 * Score for tree and forest
 **/
double treeScore(Data*, pns::node*, int);
double forestScore(Data*, std::vector<Tree*>);

/**
 * Methods to calculate split costs, or data and model costs of a
 * split separately indepenent of the type. All need two vectors
 * with the indices belonging to the two possible children.
 **/
double splitLeafData(std::vector<int>&,std::vector<int>&,Data*,pns::node*);
double splitLeafModel(std::vector<int>&,std::vector<int>&,Data*,pns::node*);
double splitLeafCosts(std::vector<int>&,std::vector<int>&,Data*,pns::node*);
double splitLeafCosts(std::vector<int>&,std::vector<int>&,Data*,pns::ad*,pns::ad*,pns::node*);
double splitLeafCosts(std::vector<int>&,std::vector<int>&,Data*,double,double,pns::node*);


/**
 * Methods to calculate the best candidate to split on and also pass the split
 * information.
 **/
pns::IndCut* splitBinary(Data*, pns::node*, int);
pns::IndCut* splitCategorial(Data*, pns::node*, int);
pns::IndCut* splitRegression(Data* d, pns::node* n, int cand);
pns::IndCut* splitNode(Data*, pns::node*, int);

/**
 * Methods to iteratively calculate costs according to a regression candidate
 * to perform in O(n log n) instead of O(n^2).
 **/
pns::IndValPair regressionScoring(Data*,pns::node*,std::vector<pns::IndValPair>&);
pns::IndValPair regressionScoringReg(Data*,pns::node*,std::vector<pns::IndValPair>&);
pns::IndValPair regressionScoringMultiway(Data*,pns::node*,std::vector<pns::IndValPair>&);
pns::IndValPair regressionScoringRegNML(Data*,pns::node*,std::vector<pns::IndValPair>&);
pns::IndValPair regressionScoringCB(Data*,pns::node*,std::vector<pns::IndValPair>&);
pns::IndValPair regressionScoringRegNML(Data* d,pns::node*,std::vector<pns::IndValPair>&);

pns::IndValPair regressionScoringInterval(Data*,pns::node*,std::vector<pns::IndValPair>&);

/**
 * Method to calculate regressor
 **/
pns::Regressor* calculateRegressor(Data*,pns::node*,std::vector<pns::IndValPair>&);
pns::Regressor* calculateQuadraticRegressor(Data*,pns::node*,std::vector<pns::IndValPair>&);

/**
 * Wrapper methods that decide which costs functions to apply
 **/
double leafData(Data*, pns::node*);
double leafModelApproximation(Data*, pns::node*);
double leafCosts(Data*, pns::node*);
double splitCategorialComplete(Data*,pns::node*,std::vector<pns::IndValPair>&);

/**
 * Generates a copy of a node with the given indices
 **/
pns::node* copyNode(pns::node*, std::vector<int>&);
pns::node* copyNode(pns::node*, std::vector<int>&, pns::Regressor*, int);

/**
 * Returns costs to choose splitPoint
 */
double getSplitPointCosts(int, int);

double binaryLeafModelApproximation(int, Data*);
double categorialLeafModelApproximationPrecal(int, int, Data*);
double categorialLeafModelApproximation(int, int, Data*);

double paramCosts(double);

#endif
