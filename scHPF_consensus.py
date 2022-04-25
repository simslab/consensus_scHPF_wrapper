#! /usr/bin/python
import argparse
import subprocess
from glob import glob
from copy import deepcopy
import joblib
import numpy as np
import igraph as ig
from sklearn import neighbors,metrics
from scipy.sparse import coo_matrix
from scipy.io import mmread
import pandas as pd
import schpf
import matplotlib as mpl
from matplotlib import pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] =42

def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-trn','--infile-train',required=True,help='Path to training data output by scHPF prep.')
    parser.add_argument('-tst','--infile-test',required=True,help='Path to test data output by scHPF prep-like.')
    parser.add_argument('-o','--outdir',required=True,help='Path to output directory.')
    parser.add_argument('-p','--prefix',required=True,help='Prefix for output files.')
    parser.add_argument('-g','--gene-infile',required=True,help='Path to gene name file.')
    parser.add_argument('-k','--k-values',required=True,type=int,nargs='+',help='Values of k to test.')
    parser.add_argument('-t','--trials',required=True,type=int,help='Number of trials for each value of k.')
    parser.add_argument('-n','--n-models',required=True,type=int,help='Number of models to consider for each value of k.')
    parser.add_argument('-m','--min-cluster-size',required=True,type=int,help='Minimum number of factors required for keeping a cluster of factors.')
    parser.add_argument('-j','--jobs',required=False,type=int,help='Maximum number of jobs to run in parallel.')
    return parser

# function for scHPF refit based on median parameters of factor clusters
def refit_local_params(X, global_params, nfactors, bp, dp, a=0.3, c=0.3, project_kw={}):
    """
    """
    project_defaults = dict(verbose=True, max_iter=50, check_freq=5)
    eta_shp, eta_rte, beta_shp, beta_rte = global_params
    
    # make a model
    eta = schpf.HPF_Gamma(np.ravel(eta_shp), np.ravel(eta_rte))
    beta = schpf.HPF_Gamma(beta_shp.T.values, beta_rte.T.values)
    model = schpf.scHPF(nfactors, eta=eta, beta=beta, bp=bp, dp=dp, a=a, c=c)           
    
    # setup projection kwarg
    for k,v in project_defaults.items():
        if k not in project_kw.keys():
            project_kw[k] = v
    loss = model.project(X, replace=True, **project_kw)
    model.loss = loss
    
    return model

# utility function for extracting model parameters into pandas
def get_param_dfs(model):
    eta_shp = pd.Series(np.ravel(model.eta.vi_shape), name=model.name)
    eta_rte = pd.Series(np.ravel(model.eta.vi_rate), name=model.name)
    beta_shp = pd.DataFrame(model.beta.vi_shape.T)
    beta_shp.index = model.name + ':' + (beta_shp.index + 1).astype(str)
    beta_rte = pd.DataFrame(model.beta.vi_rate.T, index=beta_shp.index)
    return eta_shp, eta_rte, beta_shp, beta_rte

# function for converting model parameters into pandas dataframe 
def get_spectra(models):
    eta_shp, eta_rte, beta_shp, beta_rte = zip(*[get_param_dfs(m) for m in models])
    return pd.concat(eta_shp, axis=1).T, pd.concat(eta_rte,axis=1).T, pd.concat(beta_shp), pd.concat(beta_rte)

# function for extracting gene scores from model object into pandas dataframe
def get_genescore_spectra(models):
    gene_scores = []
    for m in models:
        gs = pd.DataFrame(m.gene_score().T)
        gs.index = m.name + ':' + (gs.index + 1).astype(str)
        gene_scores.append(gs)
    return pd.concat(gene_scores)

parser = parse_user_input()
ui = parser.parse_args()

infile = ui.infile_train
outdir = ui.outdir
prefix = ui.prefix
trials = ui.trials
genes_infile = ui.gene_infile

if not ui.jobs: # if not restriction on number of parallel jobs, run an scHPF training job for each value of k
    procs=[]
    for k in ui.k_values:
        cmd = 'scHPF train -i %(infile)s -o %(outdir)s -p %(prefix)s -k %(k)d -t %(trials)d --save-all' % vars()
        p = subprocess.Popen(cmd,shell=True)
        procs.append(p)
    p_exit = [p.wait() for p  in procs]
else:  # otherwise, run only ui.jobs k-values at-a-time
    st=0
    sp=ui.jobs
    while st <= len(ui.k_values):
        procs=[]
        for k in ui.k_values[st:sp]:
            cmd = 'scHPF train -i %(infile)s -o %(outdir)s -p %(prefix)s -k %(k)d -t %(trials)d --save-all' % vars()
            p = subprocess.Popen(cmd,shell=True)
            procs.append(p)
        p_exit = [p.wait() for p in procs]
        st+=ui.jobs
        sp+=ui.jobs

# get the model objects for the top ui.n_models models for each value of k
models_str = ui.outdir+'/*scHPF*.joblib'
model_infiles = sorted(glob(models_str))
top_model_infiles = [model_infile for model_infile in model_infiles if ('reject' not in model_infile or int(model_infile.split('reject')[1][0])<ui.n_models)]
top_model_names = [model_infile.split('/')[-1].rsplit('.',-1)[0] for model_infile in top_model_infiles]
top_model_Ks = [int(model_infile.split('scHPF_K')[1].split('_')[0]) for model_infile in top_model_infiles]
top_model_dfs = pd.DataFrame(list(zip(top_model_names,top_model_Ks,top_model_infiles)),columns=['name','K','model_file'])
top_models = [joblib.load(model_infile) for model_infile in top_model_infiles]
for model,name in zip(top_models,top_model_names):
    model.name=name

for i,model in enumerate(top_models):
    if i>0:
        gscores = np.concatenate((gscores,model.gene_score()),axis=1)
    else:
        gscores=model.gene_score()

gscores = get_genescore_spectra(top_models)
gscore_cvs = (gscores.std()/gscores.mean())
top_gene_ixs = gscore_cvs.nlargest(1000).index.values
gscores = gscores[top_gene_ixs]
print(gscores.shape)

eta_shp,eta_rte,beta_shp,beta_rte = get_spectra(top_models)
eta_ex = eta_shp/eta_rte
beta_ex = beta_shp/beta_rte

# convert gene score matrix into distance matrix and knn graph
factor_dists = 1.-pd.DataFrame(gscores,index=gscores.index).T.corr().values
n_neighbors = max(5,int(0.25*len(top_model_infiles))) # heuristic for k in knn graph
adj_binary = neighbors.kneighbors_graph(factor_dists,n_neighbors,metric='precomputed')
adj = np.zeros(adj_binary.shape)
for i,j in np.stack(adj_binary.nonzero()).T:
    adj[i,j] = metrics.jaccard_score(adj_binary[i,:].A[0], adj_binary[j,:].A[0])
adj=coo_matrix(adj)

vcount = max(adj.shape)
sources,targets = adj.nonzero()
edgelist = list(zip(sources.tolist(),targets.tolist()))

# perform walktrap clustering on knn graph
knn = ig.Graph(edges=edgelist,directed=False)
knn.vs['label']=gscores.index
knn.es['width']=adj.data
knn.es['weight']=adj.data
cluster_result = knn.community_walktrap(weights=adj.data,steps=4)
nclusters = cluster_result.optimal_count+2 # heuristic for numbr of clusters
cluster_labels = pd.Series(cluster_result.as_clustering(nclusters).membership,index=gscores.index)

# compute cluster median parameters to initialization scHPF refit
eta_shp_med = eta_shp.median().values
eta_rte_med = eta_rte.median().values
beta_shp_med = beta_shp.groupby(cluster_labels).median()
beta_rte_med = beta_rte.groupby(cluster_labels).median()
eta_ex_med = eta_shp_med/eta_rte_med
beta_ex_med = beta_shp_med/beta_rte_med
min_cluster_size = ui.min_cluster_size
keep = np.where(cluster_labels.value_counts(sort=False)>=min_cluster_size)[0]
cluster_labels=cluster_labels.loc[cluster_labels.isin(keep)]

# make cluster modularity plot
pdf_outfile = outdir+'/'+prefix+'.walktrap.pdf'
optimal_count = cluster_result.optimal_count
x=np.arange(-5,10)
modularity = [cluster_result.as_clustering(optimal_count+i).modularity for i in x]
fig,ax=plt.subplots()
ax.plot(x+optimal_count,modularity)
ax.axvline(nclusters,c='r')
ax.set_xlabel('Number of Clusters')
ax.set_ylabel('Modularity')
fig.savefig(pdf_outfile)

outfile1 = outdir+'/'+prefix+'.consensus.joblib'
np.random.seed(0)
nfactors = cluster_labels.nunique()
print(nfactors)
a=0.3
c=0.3
for model in top_models:
    if model.nfactors == nfactors:
        a=model.a
        c=model.c
        break

# compute initial consensus scHPF by refitting with initialization from cluster median parameters
matrix = mmread(infile)
consensus1 = refit_local_params(matrix, (eta_shp_med, eta_rte_med,beta_shp_med.iloc[keep],beta_rte_med.iloc[keep]),
        nfactors,bp=top_models[0].bp,dp=top_models[0].dp,a=a,c=c,project_kw={'max_iter':1})
joblib.dump(consensus1,outfile1)

# compare consensus scHPF to randomly initialized model with same k using test data
outfile2 = outdir+'/'+prefix+'.consensus.final.joblib'
test_matrix = mmread(ui.infile_test)
test_loss = schpf.loss.projection_loss_function(schpf.loss.mean_negative_pois_llh,test_matrix,consensus1.nfactors,proj_kwargs={'reinit':False, 'verbose':False})

# update consensus scHPF training until convergence starting from initial consensus scHFP
consensus2 = deepcopy(consensus1)
np.random.seed(0)
consensus2.fit(matrix, loss_function=test_loss, reinit=False, verbose=True, max_iter=150)
joblib.dump(consensus2,outfile2)


