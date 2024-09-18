#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import multixrank

multixrank_obj = multixrank.Multixrank(config="dcaf.net/config_minimal.yml", wdir="dcaf.net")
# In[4]:


print(multixrank_obj)


# In[5]:


ranking_df = multixrank_obj.random_walk_rank()


# In[12]:


print(ranking_df)


# In[7]:


multixrank_obj.write_ranking(ranking_df, path="output_dcaf.net")


# In[42]:


multixrank_obj.to_sif(ranking_df, path="output_dcaf.net/dcaf.net_seedDCAF15_top3.sif", top=20)