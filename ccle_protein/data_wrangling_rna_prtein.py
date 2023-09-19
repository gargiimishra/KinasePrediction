#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd


# In[3]:


df_gene_ex = pd.read_csv("/Users/gargi/Desktop/hutch_part2/RNA_protein_analysis/CCLE RNAseq data.csv", header= 2, index_col= 1)


# In[4]:


df_gene_ex.head()


# In[5]:


df= pd.read_csv("/Users/gargi/Desktop/hutch_part2/RNA_protein_analysis/protein_quant_current_normalized.csv")


# In[6]:


df.head()


# In[ ]:


columns_to_drop = df.columns[6:49]
df = df.drop(columns=columns_to_drop)


# In[ ]:


df.head()


# In[ ]:





# In[18]:


y = df.columns


# In[19]:


y


# In[20]:


for ele in y:
    print(ele)


# In[ ]:


# find index of last '_' in the string
new_cols = []
for ele in y:
    ind = ele.rfind("_")
    ele = ele[:ind]
    new_cols.append(ele)
    


# In[ ]:


new_cols


# In[ ]:


df.columns = new_cols


# In[ ]:


df.columns


# In[8]:


df.head()


# In[ ]:




