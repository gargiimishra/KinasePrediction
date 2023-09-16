#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt


# ### Load CCLE dataset

# In[3]:


df_gene = pd.read_csv("/Users/gargi/desktop/hutch_part2/KinPredict/CCLE_RNAseq_genes_rpkm_20180929.gct/CCLE RNAseq data.csv", header= 2, index_col= 1)
df_gene.head()


# In[12]:
df_kinase= pd.read_csv("/Users/gargi/Desktop/hutch_part2/KinPredict/Data from Taran/Kinase abundance_all.csv", index_col = 0)
df_kinase.head()
gene_list = df_kinase.index
gene_list = gene_list.to_list()
print(gene_list)

filtered_ccle = df_gene.loc[gene_list]


# In[13]:


print(filtered_ccle.shape)


# ### kinase dataset_Gujral_lab

# In[8]:


df_kinase= pd.read_csv("/Users/gargi/Desktop/hutch_part2/KinPredict/Data from Taran/Kinase abundance_all.csv", index_col = 0)


# ### extracting kinase gene list from the Gujral lab dataframe

# In[30]:


gene_list = df_kinase.index
gene_list = gene_list.to_list()
gene_list, len(gene_list)


# In[81]:


df_protein= pd.read_csv("/Users/gargi/Desktop/hutch_part2/RNA_protein_analysis/protein_quant_current_normalized.csv", index_col = 1)
print(df_protein.shape)


# In[82]:


columns_to_drop = df_protein.columns[6:49]
df_protein = df_protein.drop(columns=columns_to_drop)


# In[83]:


print(df_protein.shape)


# ### renaming column names so that cell lines have same names/index

# In[84]:


y= df_protein.columns
# find index of last '_' in the string
new_cols = []
for ele in y:
    ind = ele.rfind("_")
    ele = ele[:ind]
    new_cols.append(ele)



# In[85]:


df_protein.columns = new_cols


# In[86]:


print(df_protein.columns)


# In[93]:


print(df_protein.shape)


# ### Filter the common cell lines between df_gene and df_protein

# In[21]:


df_gene_cell_lines = df_gene.columns


# In[22]:


print(df_gene_cell_lines)


# In[23]:


df_protein_cell_lines = df_protein.columns


# In[26]:


print(type(df_protein_cell_lines))


# In[25]:


print(df_protein_cell_lines)


# ### finding intersecting columns names in gene and protein

# In[94]:


common_columns = filtered_ccle.columns.intersection(df_protein.columns)
print(len(common_columns))


# In[114]:


print(common_columns)


# In[96]:


df_ccle_kinase_common_protein_cell_lines = filtered_ccle[common_columns]


# In[97]:


print(df_ccle_kinase_common_protein_cell_lines)


# In[ ]:





# In[111]:


df_protein_common_protein_cell_lines = df_protein[common_columns]


# In[123]:


print(df_protein_common_protein_cell_lines)


# In[117]:


# Remove duplicate columns pandas DataFrame
df_protein_common_protein_cell_lines = df_protein_common_protein_cell_lines.loc[:,~df_protein_common_protein_cell_lines.columns.duplicated()]


# In[129]:


df_protein_kinase_common_protein_cell_lines = df_protein_common_protein_cell_lines[df_protein_common_protein_cell_lines.index.isin(gene_list)]


# In[130]:


print(df_protein_kinase_common_protein_cell_lines)


# In[131]:


print(df_protein_kinase_common_protein_cell_lines.shape)


# In[132]:


print(df_ccle_kinase_common_protein_cell_lines.shape)


# In[ ]:
