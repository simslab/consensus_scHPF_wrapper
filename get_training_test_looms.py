#! /usr/bin/python
import argparse, loompy, sys, numpy as np
from numpy.random import choice

def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--loom', required=True,
                        help='Path to loom file.')
    parser.add_argument('-a', '--attribute', required=True,
                        help='Column attribute for sampling.')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Output prefix.')
    parser.add_argument('-n', '--n-test-cells', type=int, required=False,
                        help='Number of test cells for each sample.',
                        default=None)
    parser.add_argument('-c', '--n-test-cells-pct', type=float, required=False,
                        help='Percent of test cells for each sample.',
                        default=None)
    parser.add_argument('-s', '--seed', type=int, required=False,
                        help='Seed. Default is 0.', default=0)
    return parser

parser = parse_user_input()
ui = parser.parse_args()
np.random.seed(ui.seed)

if ui.n_test_cells is not None and ui.n_test_cells_pct is not None:
    parser.error("Options --n-test-cells and --n-test-cells-pct cannot be "
                 "supplied together. Please choose only one.")
elif ui.n_test_cells is None and ui.n_test_cells_pct is None:
    parser.error("Please specify one of the options --n-test-cells or "
                 "--n-test-cells-pct.")

try:
    data = loompy.connect(ui.loom, validate=False)
    attr = data.ca[ui.attribute]

    train_ix=[]
    test_ix=[]
    for a in set(attr):
        ix = np.where(attr == a)[0]
        if ui.n_test_cells_pct is not None and ui.n_test_cells is None:
            ix_n = len(ix)
            sample_n = int(ix_n * ui.n_test_cells_pct)
            if sample_n < 1:
                sample_n = 1
            test = choice(ix, sample_n, replace=False)
        elif ui.n_test_cells is not None and ui.n_test_cells_pct is None:
            test = choice(ix, ui.n_test_cells, replace=False)
        else:
            raise ValueError("Sampling failed. Check your --n-test-cells or"
                             "--n-test-cells-pct.")
        test_ix.extend(test)
        train_ix.extend([i for i in ix if i not in test])
    train_ix = np.array(train_ix)
    test_ix = np.array(test_ix)
    
    print(f'Original number of cells: {len(attr)}')
    print(f'Training number of cells: {len(train_ix)}')
    print(f'Testing number of cells: {len(test_ix)}')
    
    train_output = ui.prefix+'.train.loom'
    with loompy.new(train_output) as dsout:
        for (ix,selection, view) in data.scan(items=train_ix, axis=1):
            dsout.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)
    
    test_output = ui.prefix+'.test.loom'
    with loompy.new(test_output) as dsout:
        for (ix,selection, view) in data.scan(items=test_ix, axis=1):
            dsout.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)
except Exception as e:
    print("Error:", e)
    sys.exit(1)

