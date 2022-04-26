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
- scHPF (https://github.com/simslab/scHPF)

## HELPER SCRIPTS
There are a few helper scripts for sub-sampling (equal representation of cells across specimens, subjects, batches, etc) and down-sampling (equal average counts/cell across specimens, subjects, batches, etc) and gerating training and test count matrices.

For sub-sampling cells:
'''
python subsample_merged_loom.py -i INPUT_LOOM -o OUTPUT_LOOM --subsample-col-attr COL_ATTR_FOR_SUBSETTING
'''
The above will take sub-sample all of the subtypes defined by a column attribute to the same number of cells as the subtype with the fewest cells. One can sub-sample to a fixed minimum number of cells using the --number-of-cells option. One can also sub-sample a subpopulation of cells in the loom by filtering on a different column attribute using --col-attr and --col-val to identify the column attribute and the values for that attribute that you want to keep.  For example:
'''
python subsample_merged_loom.py -i INPUT_LOOM -o OUTPUT_LOOM --subsample-col-attr COL_ATTR_FOR_SUBSETTING --col-attr COL_ATTR_FOR_MASKING --col-val MASK_VALUE_1 MASK_VALUE_2
'''
This will take all of the cells with COL_ATTR_FOR_MASKING equal to either MAKSK_VALUE_1 or MASK_VALUE_2 and sub-sample those using subsets defined by COL_ATTR_FOR_SUBSETTING.

For down-sampling molecules:
'''
python downsample_loom.py -l INPUT_LOOM -a COL_ATTR_FOR_SUBSETTING -o OUTPUT_LOOM
'''
This will take all of the cells, subset them based on COL_ATTR_FOR_SUBSETTING and randomly downsample the count matrices for each subest to the same average number of molecules per cell (which will correspond to the subset with the lowest average number of molecules per cell). One can also set a fixed target average number of molecules per cell for downsampling using the --n-molecules option.

For dividing a loom file into training and test sets:
'''
python get_training_test_looms.py -l INPUT_LOOM -p PREFIX_FOR_OUTPUT_LOOMS -a COL_ATTR_FOR_SUBSETTING -n NUMBER_OF_TEST_CELLS_PER_SUBSET
'''
This will take all of the cells, subset them based on COL_ATTR_FOR_SUBSETTING, randomly select NUMBER_OF_TEST_CELLS_PER_SUBSET cells from each subset, place these cells in a loom file for the test data set, and place the remaining cells in a loom file for the training data set.

