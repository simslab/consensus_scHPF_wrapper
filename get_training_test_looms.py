#! /usr/bin/python
import loompy
import numpy as np
import argparse
from numpy.random import choice


def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l','--loom',required=True,help='Path to loom.')
    parser.add_argument('-a','--attribute',required=True,help='Column attribute for sampling.')
    parser.add_argument('-p','--prefix',required=True,help='Output prefix.')
    parser.add_argument('-n','--n-test-cells',type=int,required=True,help='Number of test cells for each sample.')
    return parser

parser = parse_user_input()
ui = parser.parse_args()

data = loompy.connect(ui.loom,validate=False)
attr = data.ca[ui.attribute]

train_ix=[]
test_ix=[]
for a in set(attr):
    ix = np.where(attr==a)[0]
    test = choice(ix,ui.n_test_cells,replace=False)
    test_ix.extend(list(test))
    train_ix.extend([i for i in ix if i not in test])
train_ix=np.array(train_ix)
test_ix = np.array(test_ix)

print(len(attr))
print(len(train_ix))
print(len(test_ix))

train_output = ui.prefix+'.train.loom'
with loompy.new(train_output) as dsout:
    for (ix,selection,view) in data.scan(items=train_ix,axis=1):
        dsout.add_columns(view.layers,col_attrs=view.ca,row_attrs=view.ra)

test_output = ui.prefix+'.test.loom'
with loompy.new(test_output) as dsout:
    for (ix,selection,view) in data.scan(items=test_ix,axis=1):
        dsout.add_columns(view.layers,col_attrs=view.ca,row_attrs=view.ra)        
