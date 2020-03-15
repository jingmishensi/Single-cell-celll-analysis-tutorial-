#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd 
import numpy as np 
import scanpy as sc


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=300)



# In[2]:


import os


# In[3]:


os.getcwd()


# In[4]:


os.chdir("C:/Users/xxj07/Desktop/R practice/python scanpy/filtered_gene_bc_matrices/")



# In[5]:


# Read  the data 
data = sc.read_10x_mtx("./hg19/", var_names='gene_symbols',                  # use gene symbols for the variable names (variables-axis index)
    cache=True)


# In[6]:


# Compute the top 25 highexpressed genes 
sc.pl.highest_expr_genes(data, n_top=25)


# In[7]:


#Filter the cells and genes according to the principle : the cells with the genes more than 200 , and the gene more than three cells are kepted 

sc.pp.filter_cells(data,min_genes=200)
sc.pp.filter_genes(data,min_cells=3)


# In[8]:


mito_genes = data.var_names.str.startswith('MT-')


# In[9]:


#Compute the mitochondria genes ratio compared to the total counts of the genes
data.obs['percent_mito'] = np.sum(
    data[:, mito_genes].X, axis=1).A1 / np.sum(data.X, axis=1).A1


# In[10]:


data.obs['n_counts'] = data.X.sum(axis=1).A1


# In[11]:


# Violin plot the n_genes , n_counts and percent_mito

sc.pl.violin(data, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True)



# In[12]:


sc.pl.scatter(data,x= "n_genes",y = "percent_mito")
sc.pl.scatter(data, x= "n_counts",y="n_genes")


# In[13]:


#Filters the cells as the following peinciple: number of genes less than 2500 ,  the percent of mitochondrial less than 0.05 

data = data[data.obs.n_genes < 2500, :]
data = data[data.obs.percent_mito < 0.05, :]


# In[14]:


sc.pp.normalize_total(data,target_sum=1e4)


# In[15]:


sc.pp.log1p(data)


# In[16]:


#Compute the high variable gene 
sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)


# In[17]:


sc.pl.highly_variable_genes(data)


# In[18]:


data = data[:, data.var.highly_variable]


# In[19]:


# Scale the data according to  the "n_counts"  and "percent_mito"
sc.pp.regress_out(data, ['n_counts', 'percent_mito'])
sc.pp.scale(data, max_value=10)


# In[20]:


sc.tl.pca(data, svd_solver='arpack')


# In[21]:


sc.pl.pca(data, color='CST3')


# In[22]:


sc.pl.pca(data)


# In[23]:


sc.pl.pca_variance_ratio(data, log=True)


# In[24]:


sc.pl.pca_variance_ratio(data)


# In[25]:


sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)


# In[26]:


sc.tl.umap(data)


# In[27]:


sc.pl.umap(data, color=['CST3', 'NKG7', 'PPBP'])


# In[28]:


sc.pl.umap(data, color=['CST3', 'NKG7', 'PPBP'], use_raw=False)


# In[29]:


sc.tl.leiden(data)


# In[30]:


sc.pl.umap(data,color = ["leiden"])


# In[31]:


sc.tl.rank_genes_groups(data, 'leiden', method='t-test')
sc.pl.rank_genes_groups(data, n_genes=25, sharey=False)


# In[ ]:



# According to the authors' github code , 
import scanpy as sc
import numpy as np

X = np.random.randint(0,1000, size= (3000,2000))
ann = sc.AnnData(np.log(X+1))
gsize = X.shape [0] / 2
ann.obs['group'] = ['a']* int (gsize) + ['b']*int (gsize)
sc.tl.rank_genes_groups(ann, 'group', method = 'wilcoxon', n_genes=2000)


# In[ ]:





# In[54]:


sc.settings.verbosity = 2


# In[497]:


sc.tl.rank_genes_groups(data, 'leiden', method='logreg')
sc.pl.rank_genes_groups(data, n_genes=25, sharey=False)


# In[32]:


new_cluster_names = [
    'CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T',
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes']
data.rename_categories('leiden', new_cluster_names)


# In[33]:


sc.pl.umap(data, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')


# In[35]:


marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']


# In[36]:


ax = sc.pl.dotplot(data, marker_genes, groupby='leiden')


# In[37]:


ax = sc.pl.stacked_violin(data, marker_genes, groupby='leiden', rotation=90)


# In[38]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc


# In[39]:


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()


# In[51]:


adata = sc.datasets.paul15()


# In[52]:


adata.X = adata.X.astype('float64')


# In[53]:


# Filtered the cells as the recipe_zheng17 
sc.pp.recipe_zheng17(adata)


# In[54]:


# PCA reduced the dimension
sc.tl.pca(adata, svd_solver='arpack')


# In[55]:


# Computing the  graph and  draw  the graph
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata)


# In[56]:


sc.pl.draw_graph(adata, color='paul15_clusters', legend_loc='on data')


# In[57]:


sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.draw_graph(adata)

#sc.tl.louvain(data)


# In[58]:


sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_diffmap')


# In[59]:


sc.tl.draw_graph(adata)


# In[61]:


sc.tl.louvain(adata)


# In[62]:


sc.tl.paga(adata, groups='louvain')


# In[63]:


sc.pl.paga(adata, color=['louvain'],title = "")


# In[64]:


sc.tl.louvain(adata, resolution=1.0)


# In[65]:



sc.tl.paga(adata, groups='louvain')


# In[66]:


import numpy as np 

sc.pl.paga(adata)


# In[67]:


sc.pl.paga(adata, color=['louvain', 'Itga2b', 'Prss34', 'Cma1'])


# In[68]:


adata.obs['louvain'].cat.categories


# In[69]:


adata.obs['louvain_anno'] = adata.obs['louvain']


# In[70]:


adata.obs['louvain_anno'].cat.categories = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10/Ery', '11', '12',
       '13', '14', '15', '16/Stem', '17', '18', '19/Neu', '20/Mk', '21', '22/Baso', '23', '24/Mo']


# In[71]:


sc.tl.paga(adata, groups='louvain_anno')


# In[72]:


sc.pl.paga(adata, threshold=0.03, show=False)


# In[73]:



sc.tl.draw_graph(adata, init_pos='paga')


# In[74]:


sc.pl.draw_graph(adata, color=['louvain_anno', 'Itga2b', 'Prss34', 'Cma1'], legend_loc='on data')


# In[75]:


pl.figure(figsize=(8, 2))
for i in range(28):
    pl.scatter(i, 1, c=sc.pl.palettes.zeileis_28[i], s=200)
pl.show()


# In[79]:


zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata.uns['louvain_anno_colors'])


# In[82]:


new_colors[[2]]


# In[83]:


new_colors[[16]] = zeileis_colors[[12]]  # Stem colors / green
new_colors[[10, 17, 5, 3, 15, 6, 18, 13, 7, 12]] = zeileis_colors[[5, 5, 5, 5, 11, 11, 10, 9, 21, 21]]  # Ery colors / red
new_colors[[20, 8]] = zeileis_colors[[17, 16]]  # Mk early Ery colors / yellow
new_colors[[4, 0]] = zeileis_colors[[2, 8]]  # lymph progenitors / grey
new_colors[[22]] = zeileis_colors[[18]]  # Baso / turquoise
new_colors[[19, 14, 2]] = zeileis_colors[[6, 6, 6]]  # Neu / light blue
new_colors[[24, 9, 1, 11]] = zeileis_colors[[0, 0, 0, 0]]  # Mo / dark blue
new_colors[[21, 23]] = zeileis_colors[[25, 25]]  # outliers / grey


# In[84]:


adata.uns['louvain_anno_colors'] = new_colors


# In[85]:


sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)


# In[86]:


adata.uns['iroot'] = np.flatnonzero(adata.obs['louvain_anno']  == '16/Stem')[0]


# In[87]:


np.flatnonzero(adata.obs['louvain_anno']  == '16/Stem')


# In[89]:


adata.obs['louvain_anno']  == '16/Stem'


# In[90]:


sc.tl.dpt(adata)


# In[91]:


gene_names = ['Gata2', 'Gata1', 'Klf1', 'Epor', 'Hba-a2',  # erythroid
              'Elane', 'Cebpe', 'Gfi1',                    # neutrophil
              'Irf8', 'Csf1r', 'Ctsg']                     # monocyte


# In[92]:


adata_raw = sc.datasets.paul15()
sc.pp.log1p(adata_raw)
sc.pp.scale(adata_raw)
adata.raw = adata_raw


# In[93]:


sc.pl.draw_graph(adata, color=['louvain_anno', 'dpt_pseudotime'], legend_loc='on data')


# In[94]:


paths = [('erythrocytes', [16, 12, 7, 13, 18, 6, 5, 10]),
         ('neutrophils', [16, 0, 4, 2, 14, 19]),
         ('monocytes', [16, 0, 4, 11, 1, 9, 24])]


# In[95]:


adata.obs['distance'] = adata.obs['dpt_pseudotime']



# In[97]:


adata.obs['clusters'] = adata.obs['louvain_anno']  # just a cosmetic change


# In[98]:


adata.uns['clusters_colors'] = adata.uns['louvain_anno_colors']


# In[ ]:




