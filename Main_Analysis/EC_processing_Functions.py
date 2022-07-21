# import
import scanpy as sc
import anndata as ad
import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as pl
import numpy as np
import scrublet as scr
import seaborn as sb
import scipy as sci
from scipy.sparse import csr_matrix
import scvelo as scv
import loompy
import leidenalg
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
from scipy import sparse
import os
import bbknn

# Single cell procession of matrix files
def scRNA(adata,resol):
    adata.var_names_make_unique()  # this is unnecessary if using 'gene_ids'
    ###Preprocessing
    #filter cells and genes
    sc.pp.filter_cells(adata, min_genes=200) 
    sc.pp.filter_genes(adata, min_cells=3)
    #filter out cells with a lot of mitochondiral DNA
    mito_genes = adata.var_names.str.startswith('MT-')
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs['pct_counts_mt'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['total_counts'] = adata.X.sum(axis=1).A1

    adata = adata[adata.obs['n_genes'] > 2000, :]
    adata = adata[adata.obs['pct_counts_mt'] < 0.1, :]

    #Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, 
    #so that counts become comparable among cells.
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=10000)
    sc.pp.log1p(adata) #logarithm
    adata.raw = adata

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    #sc.pl.highly_variable_genes(adata)

    adata = adata[:, adata.var['highly_variable']] #Actually do the filtering
    
    #Regression and PCAs

    #Scale each gene to unit variance. Clip values exceeding standard deviation 10.
    sc.pp.scale(adata, max_value=10)

    ##Principal component analysis
    #Reduce the dimensionality of the data by running principal component analysis (PCA),
    #which reveals the main axes of variation and denoises the data.
    sc.tl.pca(adata, svd_solver='arpack',n_comps=100)
    #sc.pl.pca(adata, color=['PECAM1'])

    #sc.pl.pca_variance_ratio(adata,log=True, n_pcs=100)#check how many PCAs should be considered
    
    #neighborhood graph + calculation of louvain, paga and umap
    #sc.pp.neighbors(adata, n_neighbors=10, n_pcs=77)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    ##Embedding the neighborhood graph
    #It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preservers trajectories.
    #tl for compution and pl for plotting
    sc.tl.leiden(adata, resolution=resol) #addition does clustering
    sc.tl.paga(adata,groups='leiden') #compute paga
    sc.tl.umap(adata) #compute UMAP
    #sc.tl.tsne(adata) #compute TSNE
    if "batch" in adata.obs.columns:
        sc.pl.umap(adata,color=["leiden","batch","PECAM1","PDGFRB","EPCAM"],color_map="OrRd",wspace=0.3,ncols=5)
    else:
        sc.pl.umap(adata,color=["leiden","PECAM1","PDGFRB","EPCAM"],color_map="OrRd",wspace=0.3) 
    
    import scrublet as scr
    scrub = scr.Scrublet(adata.raw.X)
    adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
    scrub.plot_histogram()
    sc.pl.umap(adata, color=["leiden", 'doublet_scores','predicted_doublets'],color_map="OrRd")
    
    return adata
# call function example
# data=scRNA(data,0.4)

# rank genes and return data frame of genes and pvalues
def Rank_Genes(adata,ident):
    sc.tl.rank_genes_groups(adata, ident, method='t-test',n_genes=300)
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
    sc.settings.verbosity = 2  # reduce the verbosity
    pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5) #Show the 10 top ranked genes per cluster
    markergenes=adata.uns['rank_genes_groups']['names']
    #Get a table with the scores and groups
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names #has a problem
    #to save in an excel file
    pdf=pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']})
    return pdf

# rank genes (specific groups) and return data frame of genes and pvalues    
def Rank_Gene_Groups(gruppe,referenz,adata,ident):#,filename):
    sc.tl.rank_genes_groups(adata, ident, groups=gruppe, reference=referenz, method='wilcoxon')
    sc.pl.rank_genes_groups(adata, groups=gruppe, n_genes=20)
    sc.pl.rank_genes_groups_violin(adata, 
                                   groups= gruppe,
                                   save='diff_expressed_genes_'+str(gruppe)+'_vs_'+str(referenz)+'.pdf',
                                   n_genes=15,jitter=False,strip=False)
    pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5) #Show the 10 top ranked genes per cluster
    markergenes=adata.uns['rank_genes_groups']['names']
    #Get a table with the scores and groups
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names #has a problem
    #to save in an excel file
    pdf=pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']})
#     pdf.to_excel('./'+filename+
#                  str(gruppe)+'_vs_'+str(referenz)+
#                  '.xlsx')
    return pdf
# call function example
# gruppe=[adata.obs.cell_types.cat.categories[5]]
# referenz=adata.obs.cell_types.cat.categories[4]
# ident='cell_types'
# df_d5_r4=Rank_Gene_Groups(gruppe,referenz,adata,ident)


# calculate cell cycle scores for umap
def cell_cycle(adata_hvg):
    cell_cycle_genes = [x.strip() for x in open('./regev_lab_cell_cycle_genes.txt')] 
    #downloaded from https://github.com/scverse/scanpy_usage/blob/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata_hvg.var_names]

    sc.tl.score_genes_cell_cycle(adata_hvg, s_genes=s_genes, g2m_genes=g2m_genes)
    adata_hvg.obs['phase'] = pd.Categorical(adata_hvg.obs['phase'])
    ref_cluster = pd.Categorical(adata_hvg.obs['phase'],categories=['S','G2M','G1','Cycling','Non-Cycling'])
    x = adata_hvg.obs['phase'] =='S'
    ref_cluster[x] = 'Cycling'
    x = adata_hvg.obs['phase'] =='G2M'
    ref_cluster[x]='Cycling'
    x = adata_hvg.obs['phase'] =='G1'
    ref_cluster[x]='Non-Cycling'
    adata_hvg.obs['proliferation'] = ref_cluster
    adata_hvg.obs['proliferation'].cat.remove_unused_categories(inplace=True)

    adata_hvg.obs['S_score']=adata_hvg.obs['S_score']
    adata_hvg.obs['G2M_score']=adata_hvg.obs['G2M_score']
    adata_hvg.uns['proliferation_colors']=['#377eb8','#bdbdbd']

    sc.pl.umap(adata_hvg,color=['proliferation','batch','phase','S_score','G2M_score'],cmap='RdGy_r',
    #           save=#'/home/srosowski/Python/write/sec_seq/day_09_12_c_m/'
    #            '_'+todayd+'_Cell_cycle_umap.svg',
              wspace=.5)
    return adata_hvg
# call function example
# adata=cell_cycle(adata)

# repeat dimensional reduction based on cell cycle regression
def cell_cycle_regression(adata):
    sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack',n_comps=100)
#     sc.pl.pca(adata, color=['PECAM1','ACTA2'])
#     sc.pl.pca_variance_ratio(adata,log=True, n_pcs=100)
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.leiden(adata,resolution=0.3) #addition does clustering
    sc.tl.paga(adata) #compute paga
    sc.tl.umap(adata)
    sc.pl.umap(adata,color=['proliferation','batch','phase','S_score','G2M_score'],cmap='RdGy_r',
#           save=#'/home/srosowski/Python/write/sec_seq/day_09_12_c_m/'
#            '_'+todayd+'_Cell_cycle_umap.svg',
          wspace=.5)
    return adata

# call function example
# adata=cell_cycle_regression(adata)


def trans_to_splice(adata, combined_loom_file_location,name_l):
    ldata = scv.read(combined_loom_file_location, cache=True)
    ldata.var_names_make_unique()
    ldata.var_names_make_unique()
    #     process ldata
    y1_name=[s[21:26] for s in ldata.obs_names]

    x=np.array(y1_name)
    for i in range(len(name_l)):
        x[np.isin(y1_name,name_l[i])]=adata.obs.batch.cat.categories[i]
    x1=[s[0:3] for s in x]
    x2=[s[3:5] for s in x]
    x=[s +''+t for s,t in zip(x1,x2)]

    y_new=[s[27:43] for s in ldata.obs_names]
    y_new_2=[s + '-1-'+t for s,t in zip(y_new,x)]
    ldata.obs_names=y_new_2
    #     copy names to loom file
    l=ldata[adata.obs_names].copy()
    #     tranfere data from adata to l
    l.obs=adata.obs
    #     l.var=adata.var
    l.uns=adata.uns
    l.obsm=adata.obsm
    #     l.varm=adata.varm
    l.obsp=adata.obsp
    return l


# calculate average expression per gene and cluster of
def  AverageGeneExpression(adata,ident):
    res = pd.DataFrame(columns=adata.raw.var_names, index=adata.obs[ident].cat.categories)                                                                                                 
    for clust in adata.obs[ident].cat.categories: 
        res.loc[clust] = adata[adata.obs[ident].isin([clust]),:].raw.X.mean(0)
    res = res.transpose()
    return res

# calculate DEGs between clusters and filter by cluster expression mean
def filtered_dotplot(adata,ident):
    sc.tl.rank_genes_groups(adata, ident, method='t-test',n_genes=300)
    sc.settings.verbosity = 2  # reduce the verbosity
    pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5) #Show the 10 top ranked genes per cluster
    markergenes=adata.uns['rank_genes_groups']['names']
    #Get a table with the groups only
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names #has a problem
    pdf=pd.DataFrame(
        {group: result[key][group]
        for group in groups for key in ['names']})
    def  AverageGeneExpression(adata,ident):
        res = pd.DataFrame(columns=adata.raw.var_names, index=adata.obs[ident].cat.categories)                                                                                                 
        for clust in adata.obs[ident].cat.categories: 
            res.loc[clust] = adata[adata.obs[ident].isin([clust]),:].raw.X.mean(0)
        res = res.transpose()
        return res
    avgedf=AverageGeneExpression(adata,ident)
    sel_deg_dict=dict()
#     sort genes by expression 
    for i in pdf.columns:
        print(i)
        x=pdf[i]
        cur_avex=avgedf[np.isin(avgedf.index,x)]#.sort_values(by=[x.name], ascending=False)
        other_col_names=list(cur_avex.columns)
        other_col_names.remove(x.name)
        cur_avex["mean_others"]=cur_avex[other_col_names].mean(axis=1)
        cur_avex["std_others"]=cur_avex[other_col_names].std(axis=1)
        cur_avex["ratio"]=cur_avex["mean_others"].div(cur_avex[x.name], axis=0)
        sel_deg=cur_avex[cur_avex.apply(lambda cu_row: abs(cu_row['ratio']) < 0.50 
                                              and cu_row[x.name] > .5 , axis=1)
                              ].sort_values(by=["ratio"], ascending=True).index

        sc.pl.dotplot(adata,
                      sel_deg,
                      groupby=ident)
        sel_deg_dict[x.name]=sel_deg
    sel_deg_df=pd.DataFrame.from_dict(sel_deg_dict, orient='index').transpose()
    return sel_deg_df


#without gene ranking
def filtered_dotplot_wo_Rank(adata,ident,pdf):
    def  AverageGeneExpression(adata,ident):
        res = pd.DataFrame(columns=adata.raw.var_names, index=adata.obs[ident].cat.categories)                                                                                                 
        for clust in adata.obs[ident].cat.categories: 
            res.loc[clust] = adata[adata.obs[ident].isin([clust]),:].raw.X.mean(0)
        res = res.transpose()
        return res
    avgedf=AverageGeneExpression(adata,ident)
    sel_deg_dict=dict()
#     sort genes by expression 
    for i in pdf.columns:
        print(i)
        x=pdf[i]
        cur_avex=avgedf[np.isin(avgedf.index,x)]#.sort_values(by=[x.name], ascending=False)
        other_col_names=list(cur_avex.columns)
        other_col_names.remove(x.name)
        cur_avex["mean_others"]=cur_avex[other_col_names].mean(axis=1)
        cur_avex["std_others"]=cur_avex[other_col_names].std(axis=1)
        cur_avex["ratio"]=cur_avex["mean_others"].div(cur_avex[x.name], axis=0)
        sel_deg=cur_avex[cur_avex.apply(lambda cu_row: abs(cu_row['ratio']) < 0.50 
                                              and cu_row[x.name] > .5 , axis=1)
                              ].sort_values(by=[x.name], ascending=False).index  #by expression in respective cluster
#                               ].sort_values(by=["ratio"], ascending=True).index # by difference to other clusters

        sc.pl.dotplot(adata,
                      sel_deg,
                      groupby=ident)
        sel_deg_dict[x.name]=sel_deg
    sel_deg_df=pd.DataFrame.from_dict(sel_deg_dict, orient='index').transpose()
    return sel_deg_df
#example
#filtered_dotplot(adata,"leiden_anno",df)


# Rank DEGs between 2 groups only
def rank_genes_groups_df(adata, group, pval_cutoff : float =None, logfc_cutoff=None): 
    d = pd.DataFrame() 
    for k in ['scores', 'names', 'logfoldchanges', 'pvals', 'pvals_adj']: 
        d[k] = adata.uns["rank_genes_groups"][k][group] 
    if pval_cutoff is not None: 
        d = d[d["pvals_adj"] < pval_cutoff] 
    if logfc_cutoff is not None: 
        d = d[d["logfoldchanges"].abs() > logfc_cutoff] 
    d['pvals_adj'] = d['pvals_adj'].apply(lambda x : 1e-320 if x == 0 else x)
    d['-log10(pvals_adj)'] = d['pvals_adj'].apply(lambda x : - np.log10(x))
    return d

# save anndata object that it can be opened in Seurat 
def write_to_csv(adata,file_name):
    raw_data = ad.AnnData(adata.raw.X)
    raw_data.obs_names=adata.raw.obs_names
    raw_data.var_names=adata.raw.var_names
    raw_data.obs=adata.obs
    raw_data.uns=adata.uns
    raw_data.obsm=adata.obsm
    raw_data.obsp=adata.obsp
    raw_data.write_csvs(file_name,
                  skip_data=False)