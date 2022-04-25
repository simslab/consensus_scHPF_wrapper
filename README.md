# consensus_scHPF_wrapper

This is set of scripts for producing a consensus scHPF model for a scRNA-seq data set in four basic steps: 0) preparing training and test data sets; 1) generation of multiple scHPF models using multiple values of k; 2) clustering of the resulting factors in 1); 3) refitting to produce a consensus scHPF model initialized from the cluster medians in 2). While steps 2) and 3) are computationally inexpensive, step 1) is rate-limiting, and so this script facilitates efficient multi-threading for step 1). 

This pipeline has the following dependencies (a few of which are somewhat inflexible with respect to version, so use of condas and/or pip is strongly recommended):

- python >= 3.7
- numba 0.50.1
- numpy 1.18.1
- scipy 1.7.0
- matplotlib
- pandas
- loompy
- scikit-learn
- igraph


