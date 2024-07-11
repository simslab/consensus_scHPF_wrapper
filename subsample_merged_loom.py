#!/usr/bin/python
import loompy
import argparse
import numpy as np

def parse_user_input():
    """
    Get and parse user input.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--output-loom',required=True,help='Path to output loom file.')
    parser.add_argument('-i','--input-loom',required=True,help='Path to input loom file.')
    parser.add_argument('-n','--number-of-cells',required=False,type=int,help='Subsample to this fixed number of cells.')
    parser.add_argument('-p','--pct-of-cells',required=False,type=float,help='Subsample this percentage of cells.')
    parser.add_argument('-c','--col-attr',required=False,help='Identify cells for sub-sampling based on this column attribute.')
    parser.add_argument('-v','--col-val',nargs='+',required=False,help='Value of column attribute for cells to send for sub-sampling.')
    parser.add_argument('-sc','--subsample-col-attr',required=True,help='Column attribute to subset cells for sub-sampling.')
    return parser

parser = parse_user_input()
ui = parser.parse_args()

np.random.seed(0)

with loompy.new(ui.output_loom) as dsout:  # Create a new, empty, loom file
    with loompy.connect(ui.input_loom,validate=False) as ds:
        print(f'Detected {ds.shape[0]} molecules across {ds.shape[1]} cells.')
        if ui.col_attr:
            cells = np.where(np.isin(ds.ca[ui.col_attr], ui.col_val))[0]
        else:
            cells = np.arange(ds.shape[1])
        
        scts = {sample: float(sum(ds.ca[ui.subsample_col_attr][cells] == sample)) for sample in set(ds.ca[ui.subsample_col_attr][cells])}
        scts_values = np.array(list(scts.values()))
        quantiles = np.quantile(scts_values, [0.25, 0.5, 0.75])
		
        print(f'Distribution of cells per {ui.subsample_col_attr}')
        print(f'count = {len(scts_values)}')
        print(f'min = {scts_values.min()}')
        print(f'mean = {scts_values.mean():.2f}')
        print(f'median = {quantiles[1]:.2f}')
        print(f'25th percentile = {quantiles[0]:.2f}')
        print(f'50th percentile (median) = {quantiles[1]:.2f}')
        print(f'75th percentile = {quantiles[2]:.2f}')
        print(f'max = {scts_values.max()}')
		
        if ui.number_of_cells:
            mnct = ui.number_of_cells
        elif ui.pct_of_cells and 0 < ui.pct_of_cells < 1:
            mnct = ui.pct_of_cells
        else:
            mnct = scts_values.min()

        if not ui.pct_of_cells:
            sfracs = {sample: mnct / scts[sample] for sample in scts}
        else:
            sfracs = {sample: mnct for sample in scts}
		
        rnd = np.random.rand(ds.shape[1])
        keep = np.array([cell for cell in cells if rnd[cell] < sfracs[ds.ca[ui.subsample_col_attr][cell]]])
        
        for (ix, selection, view) in ds.scan(items=keep, axis=1, key="Accession"):
            dsout.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)
        
        print(f'Subsampled down to {len(keep)} cells.')

