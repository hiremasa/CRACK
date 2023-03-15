#ifndef SYNTHETIC_DATA_H
#define SYNTHETIC_DATA_H

#include <vector>
#include <map>

#include "tree.h"
#include "data.h"

void generateRandomData(int, int, char, int, double, double, Data*);
double rand_normal(double, double);
int splitLeaf(Tree*, Data*, int, int, int, int);
void generateDataFromLeaves(pns::node*, Data*, int, int);

#endif
