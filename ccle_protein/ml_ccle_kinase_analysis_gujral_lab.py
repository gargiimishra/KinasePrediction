#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd

# ### Load CCLE dataset

# In[3]:


df_gene = pd.read_csv("/Users/gargi/desktop/hutch_part2/KinPredict/CCLE_RNAseq_genes_rpkm_20180929.gct/CCLE RNAseq data.csv", header= 2, index_col= 1)
df_gene.head()


# In[12]:




# In[13]:
df_kinase= pd.read_csv("/Users/gargi/Desktop/hutch_part2/KinPredict/Data from Taran/Kinase abundance_all.csv", index_col = 0)

gene_list = df_kinase.index
gene_list = gene_list.to_list()
gene_list, len(gene_list)

filtered_ccle = df_gene.loc[gene_list]

print(filtered_ccle.shape)


# ### kinase dataset_Gujral_lab

# In[8]:




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


type(df_protein_cell_lines)


# In[25]:


print(df_protein_cell_lines)


# ### finding intersecting columns names in gene and protein

# In[94]:


common_columns = filtered_ccle.columns.intersection(df_protein.columns)
len(common_columns)


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


# ### filling up null values with using .fillna(0)
#

# In[135]:


df_ccle_kinase_common_protein_cell_lines.isnull().sum()


# In[137]:


df_protein_kinase_common_protein_cell_lines.isnull().sum()


# In[139]:


df_protein_kinase_common_protein_cell_lines = df_protein_kinase_common_protein_cell_lines.fillna(0)
df_protein_kinase_common_protein_cell_lines.isnull().sum()


# ### Random forest Regressor

# In[168]:


from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error


X = df_ccle_kinase_common_protein_cell_lines.T
y = df_protein_kinase_common_protein_cell_lines.T


# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Train a random forest regressor with 100 trees
rf = RandomForestRegressor(n_estimators=100, random_state=42, max_depth=2)
rf.fit(X_train, y_train)

# Make predictions on the test set
y_pred = rf.predict(X_test)


# In[169]:


mse = mean_squared_error(y_test, y_pred)
mae = mean_absolute_error(y_test, y_pred)
print("Mean squared error: ", mse)
print("Mean absolute error:", mae)


# In[163]:


print(y_pred[3:])


# In[165]:


print(y_test[1:])


# In[170]:


from sklearn.model_selection import GridSearchCV
import numpy as np

X = df_ccle_kinase_common_protein_cell_lines.T
y = df_protein_kinase_common_protein_cell_lines.T


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Define a range of values to search for n_estimators
n_estimators = [10, 50, 100, 200]

# Number of features to consider at every split
max_features = ['auto', 'sqrt']

# Maximum number of levels in tree
#max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
max_depth=[10, 20]
max_depth.append(None)

# Minimum number of samples required to split a node
min_samples_split = [2, 5, 10]

# Minimum number of samples required at each leaf node
min_samples_leaf = [1, 2, 4]

# Method of selecting samples for training each tree
bootstrap = [True, False]

# Create the random grid
random_grid = {'n_estimators': n_estimators,
               'max_features': max_features,
               'max_depth': max_depth,
               'min_samples_split': min_samples_split,
               'min_samples_leaf': min_samples_leaf,
               'bootstrap': bootstrap}


# Create a grid search object with 5-fold cross-validation
#grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=5, scoring='neg_mean_squared_error')
grid_search = GridSearchCV(estimator=rf, param_grid=random_grid, cv=5, scoring='neg_mean_absolute_error', verbose = 2)


# Fit the grid search object to the data
grid_search.fit(X_train, y_train)

# Print the results
print("Best n_estimators:", grid_search.best_params_['n_estimators'])
print("Best MAE score:", -1 * grid_search.best_score_)


# ### Extra tree regressor

# In[ ]:


# Import the required libraries
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error

# Load your dataset
X = df_ccle_kinase_common_protein_cell_lines.T
y = df_protein_kinase_common_protein_cell_lines.T

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Create an instance of the Extra Trees Regressor class
etr = ExtraTreesRegressor(n_estimators=500, random_state=42)

# Fit the model on the training set
etr.fit(X_train, y_train)

# Predict on the test set
y_pred = etr.predict(X_test)

# Evaluate the model using mean squared error
mae = mean_absolute_error(y_test, y_pred)
print("Mean absolute error:", mae)


# In[ ]:


from sklearn.model_selection import GridSearchCV

# Define a range of values to search for n_estimators
n_estimators_range = [10, 50, 100, 200, 500]

# Load your dataset
X = df_ccle_kinase_common_protein_cell_lines.T
y = df_protein_kinase_common_protein_cell_lines.T

# Define the hyperparameters to search over
param_grid = {'n_estimators': n_estimators_range}

# Create a grid search object with 5-fold cross-validation
#grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=5, scoring='neg_mean_squared_error')
grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=5, scoring='neg_mean_absolute_error')


# Fit the grid search object to the data
grid_search.fit(X_train, y_train)

# Print the results
print("Best n_estimators:", grid_search.best_params_['n_estimators'])
print("Best MAE score:", -1 * grid_search.best_score_)

# Predict on the test set
y_pred = grid_search.predict(X_test)

# Evaluate the model using mean squared error
mae = mean_absolute_error(y_test, y_pred)
print("Mean absolute error:", mae)

