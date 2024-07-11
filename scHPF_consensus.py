#! /usr/bin/python
import argparse
import subprocess
import warnings
import joblib
import schpf
import pandas as pd
import numpy as np
import igraph as ig
import matplotlib as mpl
from sklearn import neighbors, metrics
from scipy.sparse import coo_matrix
from scipy.io import mmread
from glob import glob
from copy import deepcopy
from matplotlib import pyplot as plt

def parse_user_input():
    parser = argparse.ArgumentParser()
    parser.add_argument('-trn', '--infile-train', required=True,
                        help='Path to training data output by scHPF prep.')
    parser.add_argument('-tst', '--infile-test', required=True,
                        help='Path to test data output by scHPF prep-like.')
    parser.add_argument('-o', '--outdir', required=True,
                        help='Path to output directory.')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Prefix for output files.')
    parser.add_argument('-g', '--gene-infile', required=True,
                        help='Path to gene name file.')
    parser.add_argument('-k', '--k-values', required=True, type=int, nargs='+',
                        help='Values of k to test.')
    parser.add_argument('-t', '--trials', required=True, type=int,
                        help='Number of trials for each value of k.')
    parser.add_argument('-n', '--n-models', required=True, type=int,
                        help='Number of models to consider for each value of k.')
    parser.add_argument('-m', '--min-cluster-size', required=True, type=int,
                        help='Minimum number of factors required for keeping a cluster of factors.')
    parser.add_argument('-j', '--jobs', required=False, type=int,
                        help='Maximum number of jobs to run in parallel.')
    parser.add_argument('-r', '--re-write', required=False, default=False,
                        action='store_true',
                        help='Re-write .joblib model file. Default False.')
    parser.add_argument('-f', '--n-factors', required=False, type=int,
                        help='Number of factors for consensus scHPF. Default '
                             'is igraph VertexDendrogram.optimal_count().')
    return parser

# function for scHPF refit based on median parameters of factor clusters
def refit_local_params(X, global_params, nfactors, bp, dp, a=0.3, c=0.3,
                       project_kw={}):
    """
    """
    project_defaults = dict(verbose=True, max_iter=50, check_freq=5)
    eta_shp, eta_rte, beta_shp, beta_rte = global_params

    # make a model
    eta = schpf.HPF_Gamma(np.ravel(eta_shp), np.ravel(eta_rte))
    beta = schpf.HPF_Gamma(beta_shp.T.values, beta_rte.T.values)
    model = schpf.scHPF(nfactors, eta=eta, beta=beta, bp=bp, dp=dp, a=a, c=c)

    # setup projection kwarg
    for k, v in project_defaults.items():
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
    beta_shp.index = (model.name + '_' + str(model.nfactors) + ':' +
                      (beta_shp.index + 1).astype(str))
    beta_rte = pd.DataFrame(model.beta.vi_rate.T, index=beta_shp.index)
    return eta_shp, eta_rte, beta_shp, beta_rte

# function for converting model parameters into pandas dataframe
def get_spectra(models):
    eta_shp, eta_rte, beta_shp, beta_rte = zip(
        *[get_param_dfs(m) for m in models])
    return pd.concat(eta_shp, axis=1).T, pd.concat(eta_rte,
                                                   axis=1).T, pd.concat(
        beta_shp), pd.concat(beta_rte)

# function for extracting gene scores from model object into pandas dataframe
def get_genescore_spectra(models):
    gene_scores = []
    for m in models:
        gs = pd.DataFrame(m.gene_score().T)
        gs.index = (m.name + '_' + str(m.nfactors) + ':' +
                    (gs.index + 1).astype(str))
        gene_scores.append(gs)
    return pd.concat(gene_scores)

# parse user input
parser = parse_user_input()
ui = parser.parse_args()

infile = ui.infile_train
outdir = ui.outdir
prefix = ui.prefix
trials = ui.trials
genes_infile = ui.gene_infile
n_factors = ui.n_factors

# run scHPF training for each value of k
if not ui.jobs:  # if not restriction on number of parallel jobs, run an scHPF training job for each value of k
    procs = []
    for k in ui.k_values:
        outfile = outdir + '/' + prefix + (f'.scHPF_K{k}_*'
                                           f'{trials}trials.joblib')
        outfile = glob(outfile)
        if len(outfile) > 0 and not re_write:
            print(f'Skipping existing model for k = {k} at {outfile[0]}.')
            continue
        print(f'Running scHPF training for k = {k}')
        cmd = 'scHPF train -i %(infile)s -o %(outdir)s -p %(prefix)s -k %(k)d -t %(trials)d --save-all' % vars()
        p = subprocess.Popen(cmd, shell=True)
        procs.append(p)
    p_exit = [p.wait() for p in procs]
else:  # otherwise, run only ui.jobs k-values at-a-time
    st = 0
    sp = ui.jobs
    while st <= len(ui.k_values):
        procs = []
        for k in ui.k_values[st:sp]:
            outfile = outdir + '/' + prefix + (f'.scHPF_K{k}_*'
                                           f'{trials}trials.joblib')
            outfile = glob(outfile)
            if len(outfile) > 0 and not ui.re_write:
                print(f'Skipping existing model for k = {k} at {outfile[0]}.')
                continue
            print(f'Running scHPF training for k = {k}')
            cmd = 'scHPF train -i %(infile)s -o %(outdir)s -p %(prefix)s -k %(k)d -t %(trials)d --save-all' % vars()
            p = subprocess.Popen(cmd, shell=True)
            procs.append(p)
        p_exit = [p.wait() for p in procs]
        st += ui.jobs
        sp += ui.jobs

# get the model objects for the top ui.n_models models for each value of k
models_str = ui.outdir + '/*scHPF*.joblib'
model_infiles = sorted(glob(models_str))
top_model_infiles = [model_infile for model_infile in model_infiles if (
        'reject' not in model_infile or int(model_infile.split('reject')[1][
                                                0]) < ui.n_models)]
top_model_names = []
for model_infile in top_model_infiles:
    model_name = model_infile.split('/')[-1].rsplit('.', -1)[0]
    if 'reject' in model_infile:
        rejectn = model_infile.split('/')[-1].split('.')[1].split('_')[-1]
        top_model_names.append(f'{model_name}_{rejectn}')
    else:
        top_model_names.append(model_name)

top_model_Ks = [int(model_infile.split('scHPF_K')[1].split('_')[0]) for
                model_infile in top_model_infiles]
top_model_dfs = pd.DataFrame(list(zip(top_model_names, top_model_Ks,
                                      top_model_infiles)),
                             columns=['name', 'K', 'model_file'])
top_models = [joblib.load(model_infile) for model_infile in top_model_infiles]
for model, name in zip(top_models, top_model_names):
    model.name = name

# concatenate factors across models and select highly variable genes
gscores = get_genescore_spectra(top_models)
print(f'\nClustering {gscores.shape[0]} factors across {len(top_models)} '
      f'models with {gscores.shape[1]} genes')
gscore_cvs = (gscores.std() / gscores.mean())
print(f'Gene coefficient of variation distribution\n'
      f'{gscore_cvs.describe()}\n'
      f'Selecting top 1000 HVGs for clustering')
top_gene_ixs = gscore_cvs.nlargest(1000).index.values
gscores = gscores[top_gene_ixs]

# convert gene score matrix into knn graph
n_neighbors = max(5, int(0.25 * len(
    top_model_infiles)))  # heuristic for k in knn graph
adj_binary = neighbors.kneighbors_graph(pd.DataFrame(gscores,
                                                     index=gscores.index),
                                        n_neighbors, metric='euclidean')
adj = np.zeros(adj_binary.shape)
for i, j in np.stack(adj_binary.nonzero()).T:
    adj[i, j] = metrics.jaccard_score(adj_binary[i, :].A[0],
                                      adj_binary[j, :].A[0])

adj = coo_matrix(adj)
sources, targets = adj.nonzero()
edgelist = list(zip(sources.tolist(), targets.tolist()))

# perform walktrap clustering on knn graph
knn = ig.Graph(edges=edgelist, directed=False)
knn.vs['label'] = gscores.index
knn.es['width'] = adj.data
knn.es['weight'] = adj.data
cluster_result = knn.community_walktrap(weights=adj.data, steps=4)

# Select number of clusters based on modularity
max_k = len(cluster_result.merges)
min_k = 1 # find the minimum valid number of clusters
for i in range(1, max_k + 1):
    try:
        cluster_result.as_clustering(i)
        min_k = i
        break
    except ig._igraph.InternalError:
        continue

if n_factors is None: # Check for user-supplied number of factors
    nclusters = cluster_result.optimal_count
elif 0 < n_factors <= max_k:
    nclusters = n_factors
else:
    warnings.warn(f'"--n-factors" must be between 1 and {max_k}. '
                  f'Using default value.')
    nclusters = cluster_result.optimal_count

if nclusters < min_k:
    warnings.warn(f'Optimal number of clusters is {nclusters}, but minimum '
                  f'number of clusters where modularity can be calculated is '
                  f'{min_k}. Using minimum valid number of clusters.')
    nclusters = min_k

print(f'Number of clusters: {nclusters}')
cluster_labels = pd.Series(cluster_result.as_clustering(nclusters).membership,
                           index=gscores.index)
min_cluster_size = ui.min_cluster_size
keep = np.where(cluster_labels.value_counts(sort=False) >= min_cluster_size)[0]
cluster_labels = cluster_labels.loc[cluster_labels.isin(keep)]
print(f'Number of clusters with n_factors > {min_cluster_size}: {nclusters}')

# compute modularity of Walktrap clustering across all posible clusters
ks = np.arange(min_k, max_k + 1)
print(f'\nCalculating modularity for {min_k} - {max_k + 1} k clusters')
modularity = []
for i in ks:
    try:
        modularity.append(cluster_result.as_clustering(i).modularity)
    except ig._igraph.InternalError as e:
        print(f"Error calculating modularity at k {i}: {e}")

modularity = pd.DataFrame({'k': ks, 'modularity': modularity})
modularity['diff_modularity'] = modularity['modularity'].diff().fillna(0)
modularity.to_csv(outdir + '/' + prefix + '.walktrap.tsv', sep='\t',
                  index=False)

# plot cluster modularity
pdf_outfile = outdir + '/' + prefix + '.walktrap.pdf'
step = min(len(modularity)//10, len(modularity)//3)
plt.figure(figsize=(9, 6))
plt.rcParams.update({'font.size': 14})
plt.plot(modularity['k'], modularity['modularity'])
plt.xlim(modularity['k'].min(), modularity['k'].max())
plt.xticks(modularity['k'].iloc[::step])
plt.axvline(nclusters, c='r')
plt.text(nclusters + 0.5, plt.ylim()[1] * 0.6, f'Optimal k = {nclusters}',
         color='r')
plt.title(f'Modularity of Walktrap across {gscores.shape[0]}\n'
          f'factors from n_models ({top_model_dfs["name"].nunique()}) * '
          f'k ({top_model_dfs["K"].nunique()}) models')
plt.xlabel('Number of Clusters (k)')
plt.ylabel('Modularity')
plt.grid(True)
plt.tight_layout()
plt.savefig(pdf_outfile)
plt.close()

# compute cluster median parameters to initialize scHPF refit
eta_shp, eta_rte, beta_shp, beta_rte = get_spectra(top_models)
eta_ex = eta_shp / eta_rte
beta_ex = beta_shp / beta_rte
eta_shp_med = eta_shp.median().values
eta_rte_med = eta_rte.median().values
beta_shp_med = beta_shp.groupby(cluster_labels).median()
beta_rte_med = beta_rte.groupby(cluster_labels).median()
eta_ex_med = eta_shp_med / eta_rte_med
beta_ex_med = beta_shp_med / beta_rte_med

# initializing consensus scHPF by refitting with cluster median parameters
outfile1 = outdir + '/' + prefix + '.consensus.joblib'
np.random.seed(0)
nfactors = cluster_labels.nunique()
assert nfactors > 0, 'No valid clusters found.'
print(f'\nComputing consensus scHPF with {nfactors} factors')
a = 0.3
c = 0.3
for model in top_models:
    if model.nfactors == nfactors:
        a = model.a
        c = model.c
        break

matrix = mmread(infile)
consensus1 = refit_local_params(matrix, (eta_shp_med, eta_rte_med,
                                         beta_shp_med.iloc[keep],
                                         beta_rte_med.iloc[keep]), nfactors,
                                bp=top_models[0].bp, dp=top_models[0].dp,
                                a=a, c=c, project_kw={'max_iter': 1})
joblib.dump(consensus1, outfile1)

# compare consensus scHPF to randomly initialized model with same k using test data
outfile2 = outdir + '/' + prefix + '.consensus.final.joblib'
test_matrix = mmread(ui.infile_test)
test_loss = schpf.loss.projection_loss_function(
    schpf.loss.mean_negative_pois_llh, test_matrix, consensus1.nfactors,
    proj_kwargs={'reinit': False, 'verbose': False})

# update consensus scHPF training until convergence starting from initial consensus scHFP
consensus2 = deepcopy(consensus1)
np.random.seed(0)
consensus2.fit(matrix, loss_function=test_loss, reinit=False, verbose=True,
               max_iter=150)
joblib.dump(consensus2, outfile2)
print('Consensus scHPF has finished.')

