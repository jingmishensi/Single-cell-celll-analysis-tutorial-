#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
from matplotlib import rcParams


# In[3]:


sc.set_figure_params(dpi=80, color_map='viridis')
sc.settings.verbosity = 2
sc.logging.print_versions()


# In[4]:


pbmc = sc.datasets.pbmc68k_reduced()


# In[6]:


pbmc.var


# In[7]:


pbmc.obs


# In[8]:


pbmc.uns


# In[9]:


pbmc


# In[11]:


rcParams['figure.figsize'] = 4, 4


# In[12]:


sc.pl.umap(pbmc, color=['bulk_labels'], s=50)


# In[13]:


marker_genes = ['CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ',  'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
                'FCGR3A', 'FCER1A', 'CST3']


# In[14]:


ax = sc.pl.stacked_violin(pbmc, marker_genes, groupby='bulk_labels')
                


# In[15]:


ax = sc.pl.stacked_violin(pbmc, marker_genes)


# In[17]:


ax = sc.pl.stacked_violin(pbmc, marker_genes, groupby='bulk_labels',
                         var_group_positions=[(7, 8)], var_group_labels=['NK'])


# In[20]:


ax = sc.pl.stacked_violin(pbmc, marker_genes, groupby='bulk_labels', swap_axes=True,
                         var_group_positions=[(7, 8)], var_group_labels=['NK'], dendrogram=True)


# In[21]:


marker_genes_dict = {'B-cell': ['CD79A', 'MS4A1'],
                     'T-cell': 'CD3D',
                     'T-cell CD8+': ['CD8A', 'CD8B'],
                     'NK': ['GNLY', 'NKG7'],
                     'Myeloid': ['CST3', 'LYZ'],
                     'Monocytes': ['FCGR3A'],
                     'Dendritic': ['FCER1A']}


# In[22]:


ax = sc.pl.dotplot(pbmc, marker_genes_dict)


# In[23]:


ax = sc.pl.dotplot(pbmc, marker_genes_dict, groupby='bulk_labels')


# In[24]:


ax = sc.pl.dotplot(pbmc, marker_genes, groupby='bulk_labels', dendrogram=True, dot_max=0.5, dot_min=0.3, standard_scale='var')


# In[28]:


ax = sc.pl.dotplot(pbmc, marker_genes_dict, groupby='bulk_labels', dendrogram=True,
                   standard_scale='var', smallest_dot=40, color_map='Blues', figsize=(8,5))


# In[29]:


ax = sc.pl.dotplot(pbmc, marker_genes, groupby='louvain',
              var_group_positions=[(0,1), (11, 12)],
              var_group_labels=['B cells', 'dendritic'],
              figsize=(12,4), var_group_rotation=0, dendrogram='dendrogram_louvain')


# In[33]:


ax = sc.pl.dotplot(pbmc, marker_genes, groupby='louvain',
              var_group_positions=[(0,1), (11, 12)],
              var_group_labels=['B cells', 'dendritic'],
              figsize=(12,4))


# In[34]:


pbmc.uns


# In[36]:


ax = sc.pl.dotplot(pbmc, marker_genes, groupby='louvain',
              var_group_positions=[(0,1), (11, 12)],
              var_group_labels=['B cells', 'dendritic'],
              figsize=(12,4), var_group_rotation=0, dendrogram='dendrogram_louvain')


# In[37]:


ax = sc.pl.dotplot(pbmc, marker_genes, groupby='louvain',
              var_group_positions=[(0,1), (11, 12)],
              var_group_labels=['B cells', 'dendritic'],
              figsize=(12,4), var_group_rotation=1, dendrogram='dendrogram_louvain')


# In[38]:


ax = sc.pl.dotplot(pbmc, marker_genes, groupby='louvain',
              var_group_positions=[(0,1), (11, 12)],
              var_group_labels=['B cells', 'dendritic'],
              figsize=(12,4), var_group_rotation=0, dendrogram='dendrogram_louvain')


# In[39]:


pbmc.uns


# In[40]:


gs = sc.pl.matrixplot(pbmc, marker_genes_dict, groupby='bulk_labels')


# In[41]:


gs = sc.pl.matrixplot(pbmc, marker_genes_dict, groupby='bulk_labels', dendrogram=True, standard_scale='var')


# In[ ]:


gs = sc.pl.matrixplot(pbmc, marker_genes_dict, groupby='bulk_labels', dendrogram=True, standard_scale='var')


# In[42]:


help(sc.pl.matrix)


# In[43]:


sc.pl.matrix


# In[44]:


help(sc.pl.dotplot)


# In[46]:


marker_genes_2 = [x for x in marker_genes if x in pbmc.var_names]


# In[47]:


marker_genes_2


# In[48]:


gs = sc.pl.matrixplot(pbmc, marker_genes_2, groupby='bulk_labels', dendrogram=True,
                      use_raw=False, vmin=-3, vmax=3, cmap='bwr',  swap_axes=True, figsize=(5,6))


# In[50]:


# swap_axes=True means exchange the x and y axis
gs = sc.pl.matrixplot(pbmc, marker_genes_2, groupby='bulk_labels', dendrogram=True,
                      use_raw=False, vmin=-3, vmax=3, figsize=(5,6))


# In[51]:


ax = sc.pl.heatmap(pbmc,marker_genes_dict, groupby='louvain')


# In[52]:


ax = sc.pl.heatmap(pbmc, marker_genes, groupby='louvain', figsize=(5, 8),
              var_group_positions=[(0,1), (11, 12)], use_raw=False, vmin=-3, vmax=3, cmap='bwr',
              var_group_labels=['B cells', 'dendritic'], var_group_rotation=0, dendrogram='dendrogram_louvain')


# In[53]:


import numpy as np
ad = pbmc.copy()
ad.raw.X.data = np.exp(ad.raw.X.data)


# In[54]:


ax = sc.pl.tracksplot(ad,marker_genes, groupby='louvain')


# In[56]:



sc.pl.rank_genes_groups_dotplot(pbmc, n_genes=4)


# In[57]:


sc.pl.rank_genes_groups_dotplot(pbmc, n_genes=4)


# In[58]:


axs = sc.pl.rank_genes_groups_dotplot(pbmc, groupby='louvain', n_genes=4, dendrogram='dendrogram_louvain')


# In[60]:


axs = sc.pl.rank_genes_groups_matrixplot(pbmc, n_genes=3, standard_scale='var', cmap='Blues')


# In[61]:


axs = sc.pl.rank_genes_groups_matrixplot(pbmc, n_genes=3, use_raw=False, vmin=-3, vmax=3, cmap='bwr')


# In[64]:


sc.pl.rank_genes_groups_stacked_violin(ad, n_genes=3)


# In[65]:


# setting row_palette='slateblue' makes all violin plots of the same color
sc.pl.rank_genes_groups_stacked_violin(ad, n_genes=3, row_palette='slateblue')


# In[66]:


sc.pl.rank_genes_groups_stacked_violin(ad, n_genes=3, swap_axes=True, figsize=(6, 10), width=0.4)


# In[67]:


sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=3, standard_scale='var')


# In[68]:


sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=3)


# In[69]:


sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=3, use_raw=False, swap_axes=True, vmin=-3, vmax=3, cmap='bwr')


# In[70]:



sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=10, use_raw=False, swap_axes=True, show_gene_labels=False,
                                vmin=-3, vmax=3, cmap='bwr')


# In[71]:



sc.pl.rank_genes_groups_heatmap(pbmc, n_genes=10, use_raw=False, swap_axes=True,
                                vmin=-3, vmax=3, cmap='bwr')


# In[72]:


sc.pl.rank_genes_groups_tracksplot(ad, n_genes=3)


# In[73]:


sc.pl.rank_genes_groups_violin(pbmc,  n_genes=5, jitter=False)


# In[74]:


sc.pl.rank_genes_groups_violin(pbmc)


# In[75]:


ax = sc.pl.dendrogram(pbmc, 'bulk_labels')


# In[76]:


sc.tl.dendrogram(pbmc, 'bulk_labels', var_names=marker_genes, use_raw=True)


# In[77]:


pbmc.raw


# In[78]:


ax = sc.pl.dendrogram(pbmc, 'bulk_labels', orientation='left')


# In[79]:


sc.tl.dendrogram(pbmc, 'bulk_labels', n_pcs=30)


# In[80]:


ax = sc.pl.correlation_matrix(pbmc, 'bulk_labels')


# In[ ]:




