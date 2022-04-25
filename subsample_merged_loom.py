#! /usr/bin/python
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
	parser.add_argument('-c','--col-attr',required=False,help='Identify cells for sub-sampling based on this column attribute.')
	parser.add_argument('-v','--col-val',nargs='+',required=False,help='Value of column attribute for cells to send for sub-sampling.')
	parser.add_argument('-sc','--subsample-col-attr',required=True,help='Column attribute to subset cells for sub-sampling.')
	return parser

parser = parse_user_input()
ui = parser.parse_args()

with loompy.new(ui.output_loom) as dsout:  # Create a new, empty, loom file
	with loompy.connect(ui.input_loom) as ds:
		if ui.col_attr:
			cells = np.where(np.isin(ds.ca[ui.col_attr],ui.col_val))[0]
		else:
			cells = np.arange(ds.shape[1])
		scts = {sample:float(sum(ds.ca[ui.subsample_col_attr][cells]==sample)) for sample in set(ds.ca[ui.subsample_col_attr][cells])}
		if not ui.number_of_cells:
			mnct = np.min(list(scts.values()))
		else:
			mnct = ui.number_of_cells
		sfracs = {sample:mnct/scts[sample] for sample in set(ds.ca[ui.subsample_col_attr][cells])}
		rnd = np.random.rand(ds.shape[1])
		keep = np.array([cell for cell in cells if rnd[cell] < sfracs[ds.ca[ui.subsample_col_attr][cell]]]) 
		for (ix, selection, view) in ds.scan(items=keep, axis=1, key="Accession"):
			dsout.add_columns(view.layers,col_attrs=view.ca,row_attrs=view.ra)
