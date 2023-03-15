### Copyright (c) 2018  Alexander Marx  [amarx@mpi-inf.mpg.de]
### All rights reserved.  See the file COPYING for license terms. 

## Background information

The Crack algorithm was published in the following paper:

**Causal Inference on Multivariate and Mixed-Type Data**

Given data over the joint distribution of two random variables X and Y , we consider the problem of inferring the most likely causal direction between
X and Y. In particular, we consider the general case where both X and Y may be univariate or multivariate, and of the same or mixed data types. We
take an information theoretic approach, based on Kolmogorov complexity, from which it follows that first describing the data over cause and then that
of effect given cause is shorter than the reverse direction.
The ideal score is not computable, but can be approximated through the Minimum Description Length (MDL) principle. Based on MDL, we propose
two scores, one for when both X and Y are of the same single data type, and one for when they are mixed-type. We model dependencies between
X and Y using classification and regression trees. As inferring the optimal model is NP-hard, we propose Crack, a fast greedy algorithm to determine
the most likely causal direction directly from the data.
Empirical evaluation on a wide range of data shows that Crack reliably, and with high accuracy, infers the correct causal direction on both univariate
and multivariate cause- effect pairs over both single and mixed-type data.

If you use the Crack algorithm please cite the according paper.

## Citation

@misc{marx:18:crack,
    author = {Alexander Marx and Jilles Vreeken},
    title = {{Causal Inference on Multivariate and Mixed-Type Data}},
    booktitle={Joint European Conference on Machine Learning and Knowledge Discovery in Databases},
    year={2018},
    organization={Springer}
}

# Installation & execution
Enter the "code" folder and run
make
make cp     //copies executable into testPackage folder

To plot output trees graphviz is required.

To run the code write

crack -h

for an explanation of the arguments that can be passed.

**Run Tübingen univariate pairs:**
Go to the folder "testPackage" and run

bash run_tuebingen.sh

To change the location of the output and input please edit the script.
The results table is ordered the followin:
NameOfTestSet   Confidence  CausalDirection Runtime

**Run multivariate pais:**
Go to folder "testPackage" and run

bash run_mv.sh

In the script output and input parametes can be changed.
The results table is ordered the followin:
NameOfTestSet   Confidence  CausalDirection Runtime

**Synthetic data**
To generate synthetic data, we implemented an R script called "test_synthetic_data.R". The script
"apply_synthetic_tests.R" should be executed step by step to obtain the results for synthetic data
and obtain the decision rate plot. All these secript are located in the folder "testPackage".

