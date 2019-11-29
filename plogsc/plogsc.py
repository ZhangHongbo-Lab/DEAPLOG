from __future__ import division
from scipy.stats import fisher_exact, hypergeom,binom
from random import shuffle
import operator
import random

import math
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from scipy import sparse
from sympy import *
import matplotlib.pyplot as plt

from multiprocessing import Pool
from functools import partial

import datetime  #change 1
import fisher

def get_DEG_single(rdata:'raw annoData',
                   adata:'annoData',
                   group_key:'cell type of annoData',
                   ratio = 0.2,
                   p_threshold=0.001,
                   q_threshold=0.01):
    
    """
    get the differentially expressed genes of each cell cluster. these genes is uniq in each cell cluster.
    
    Parameters
    ----------
    rdata : raw annoData,we recommoned you use the normalized raw data
    adata : annoData, the annoData filtered by HVG is recommended;
    group_key : the cluster of cell,default,'louvain',also 'leiden' or cell_type which defined by user;
    ratio : the possibility whether gene is a  differentially expressed gene. the value is between 0 and 1.
    p_threshold : the threshold of p-value. the value is between 0 and 1.
    q_threshold : the threshold of q_value. the value is between 0 and 1.
    
    Returns
    -------
    a dataFrame of genes which is differentially expressed in each cell type and uniq in each cell type,
    which contains cell type,gene name ,ratio, p-value and q-value.  
    """
    
    #start = datetime.datetime.now()
    #print(start);
    
    
    # get the raw data frame
    print('get the raw data frame...')
    rdata_df = rdata.to_df();
    num_allCells = len(rdata_df.index)
    
    # struct the cell type sets for enrichment analysis
    print('struct the cell type sets for enrichment analysis...')
    adata_cell_type_df = pd.DataFrame(adata.obs[group_key]);
    cell_type_index = pd.Categorical(adata.obs[group_key]).categories;
    
    cell_sets = dict();
    
    for ct in cell_type_index:
        cell_sets[ct] = list(adata_cell_type_df.loc[adata_cell_type_df[group_key] ==ct,:].index);
    
    iter_genes = rdata_df.columns;
    
    print('Fisher_test_for_each_gene...')
      
    genes_ratio = dict();
    genes_pv = dict();
    genes_qv =dict();
    num = 0;
    for gene in iter_genes:
        num +=1; 
        gene_ratio, gene_pv, gene_qv = Fisher_test_for_each_gene(rdata_df,cell_sets,num_allCells,gene);
        if gene_ratio[gene] == []:
            continue;
        else:
            genes_ratio = dict(genes_ratio,**gene_ratio);
            genes_pv = dict(genes_pv,**gene_pv);
            genes_qv = dict(genes_qv,**gene_qv);
        if num%100 ==0:
            print('whole ',num,' genes have been done.')

    print('merge differentially expressed genes...')
    genes_ratio_df = pd.DataFrame(genes_ratio);
    genes_ratio_df.index = cell_type_index;
    
    genes_pv_df = pd.DataFrame(genes_pv);
    genes_pv_df.index = cell_type_index;
    
    genes_qv_df = pd.DataFrame(genes_qv);
    genes_qv_df.index = cell_type_index;
    
    ra = ratio;
    pt = p_threshold;
    qt = q_threshold;

    ct_gene_r_p_q_list =[];
    
    for gene in genes_ratio_df.columns:
        gene = gene;
        gene_on_ct = genes_ratio_df.loc[:,gene].idxmax(axis=0)
        
        if genes_ratio_df.loc[gene_on_ct,gene]>=ra:
            
            if  genes_pv_df.loc[gene_on_ct,gene]<=pt and genes_qv_df.loc[gene_on_ct,gene]<=qt:
                ct = gene_on_ct;
                g = gene;
                r = genes_ratio_df.loc[gene_on_ct,gene];
                p = genes_pv_df.loc[gene_on_ct,gene];
                q = genes_qv_df.loc[gene_on_ct,gene];
                
                ct_gene_r_p_q_list.append([ct,g,r,p,q]);
            else:
                continue;
        else:
            continue;
    markers_s = pd.DataFrame(ct_gene_r_p_q_list,columns=['cell_type','gene_name','ratio','p_value','q_value']);
    print('Done!')
    return markers_s

def get_DEG_multiple(rdata:'raw annoData',
                     adata:'annoData',
                     group_key:'cell type of annoData',
                     ratio=0.2,
                     p_threshold=0.001,
                     q_threshold=0.01):
    
   
    """
    Get the differentially expressed genes of each cell cluster. But these genes is not uniq in each cell cluster, 
    which means a gene could be differentially expressed in two or more cell clusters.
    
    Parameters
    ----------
    rdata : raw annoData, we recommoned you use the normalized raw data;
    adata : annoData, the annoData filtered by HVG is recommended;
    group_key : the cluster of cell,default,'louvain',also 'leiden' or cell_type which defined by user;
    ratio : the possibility whether gene is a  differentially expressed gene. the value is between 0 and 1;
    p_threshold : the threshold of p-value. the value is between 0 and 1;
    q_threshold : the threshold of q_value. the value is between 0 and 1;
    
    Returns
    -------
    a dataFrame of cell_clusters which contain differentially expressed genes, keys are the cell types and 
    values are arrays which contain gene name ,ratio, p-value and q-value.   
    """
    
    #start = datetime.datetime.now()
    #print(start);
    
    
    
    # get the raw data frame
    print('get the raw data frame...')
    rdata_df = rdata.to_df();
    num_allCells = len(rdata.obs_names)
    
    # struct the cell type sets for enrichment analysis
    print('struct the cell type sets for enrichment analysis...')
    adata_cell_type_df = pd.DataFrame(adata.obs[group_key]);
    cell_type_index = pd.Categorical(adata.obs[group_key]).categories;
    
    cell_sets = dict();
    
    for ct in cell_type_index:
        cell_sets[ct] = list(adata_cell_type_df.loc[adata_cell_type_df[group_key] ==ct,:].index);
    
    iter_genes = rdata_df.columns;
    
    print('Fisher_test_for_each_gene...')
      
    genes_ratio = dict();
    genes_pv = dict();
    genes_qv =dict();
    num = 0;
    for gene in iter_genes:
        num +=1; 
        gene_ratio, gene_pv, gene_qv = Fisher_test_for_each_gene(rdata_df,cell_sets,num_allCells,gene);
        if gene_ratio[gene] == []:
            continue;
        else:
            genes_ratio = dict(genes_ratio,**gene_ratio);
            genes_pv = dict(genes_pv,**gene_pv);
            genes_qv = dict(genes_qv,**gene_qv);
        if num%100 ==0:
            print('whole ',num,' genes have been done.')

    print('merge differentially expressed genes...')
    genes_ratio_df = pd.DataFrame(genes_ratio);
    genes_ratio_df.index = cell_type_index;
    genes_ratio_df = genes_ratio_df.T;
    
    genes_pv_df = pd.DataFrame(genes_pv);
    genes_pv_df.index = cell_type_index;
    genes_pv_df = genes_pv_df.T;
    
    genes_qv_df = pd.DataFrame(genes_qv);
    genes_qv_df.index = cell_type_index;
    genes_qv_df = genes_qv_df.T;
    
    ra = ratio;
    pt = p_threshold;
    qt = q_threshold;
    deg_genes_ratio_p_q_list = [];

    for ct in cell_type_index:
        up_ratio_ct_genes_list = list(genes_ratio_df.loc[genes_ratio_df[ct]>=ra,:].index);
        down_pv_ct_genes_list = list(genes_pv_df.loc[genes_pv_df[ct]<=pt,:].index);
        down_qv_ct_genes_list = list(genes_qv_df.loc[genes_qv_df[ct]<=qt,:].index);
        
        ct_deg_genes_list = list(set(up_ratio_ct_genes_list)&set(down_pv_ct_genes_list)&set(down_qv_ct_genes_list));
        
        if ct_deg_genes_list ==[]:
            continue;
            
        else:
            for gene in ct_deg_genes_list:
                ct=ct;
                g = gene;
                r = genes_ratio_df.loc[gene,ct];
                p = genes_pv_df.loc[gene,ct];
                q = genes_qv_df.loc[gene,ct];
                deg_genes_ratio_p_q_list.append([ct,g,r,p,q]);
    
    markers_m = pd.DataFrame(deg_genes_ratio_p_q_list,columns=['cell_type','gene_name','ratio','p_value','q_value']);
    print('Done!')
    return markers_m

def get_genes_location_pseudotime(rdata:'annoData of rdata',
                                  adata: 'annoData of adata',
                                  group_key:'cell type of annoData',
                                  gene_matrix:'fileName of markers matrix from results by get_DEG_single or get_DEG_multiple',
                                  obsm:'obsm'):
    
    """
    Get the differentially expressed genes of each cell cluster. But these genes is not uniq in each cell cluster, 
    which means a gene could be differentially expressed in two or more cell clusters.
    
    Parameters
    ----------
    rdata : raw annoData, we recommoned you use the normalized raw data;
    adata : annoData, the annoData filtered by HVG is recommended;
    group_key : the cluster of cell,default,'louvain',also 'leiden' or cell_type which defined by user;
    gene_matrix:'fileName of markers matrix from results by get_DEG_single or get_DEG_multiple';
    obsm:'obsm',
    n_process: 'int';
    
    Returns
    -------
    a dict of cell_clusters which contain differentially expressed genes, keys are the cell types and 
    values are arrays which contain gene name ,ratio, p-value and q-value.   
    """    
    start = datetime.datetime.now()
    #print(start);
    
    #pool = Pool(n_process);
    
    rdata_df = rdata.to_df();
    
    adata_obsm_df = pd.DataFrame(adata.obsm[obsm]);
    adata_obsm_df.index = adata.obs_names;
    adata_obsm_df['dpt_pseudotime'] = adata.obs.dpt_pseudotime;
    
    markers_s_LGPS = pd.DataFrame(adata.uns[gene_matrix],
                              columns=['cell_type','gene_name','ratio','p_value','q_value']);
    markers_s_LGPS[['ratio','p_value','q_value']] = markers_s_LGPS[['ratio','p_value','q_value']].apply(pd.to_numeric);
    markers_s_LGPS_rdata = markers_s_LGPS.loc[markers_s_LGPS['gene_name'].isin(list(rdata.var_names)),];
    #markers_s_LGPS_rdata.index = list(markers_s_LGPS_rdata['gene_name']);
    
    adata_cell_type_df = pd.DataFrame(adata.obs[group_key]);
    cell_type_index = pd.Categorical(adata.obs[group_key]).categories;
    
    cell_sets = dict();
    
    for ct in cell_type_index:
        cell_sets[ct] = list(adata_cell_type_df.loc[adata_cell_type_df[group_key] ==ct,:].index);

    iter_genes = rdata_df.columns;
    gene_pseudotime_locates =list();
    num = 0;
    for g1 in  iter_genes:
        num +=1; 
        g1_pseudotime_locates = calculate_genes_pseudotime_location(rdata_df,
                                                                    cell_sets,
                                                                    markers_s_LGPS_rdata,
                                                                    adata_obsm_df,
                                                                    g1);
        if g1_pseudotime_locates == []:
            continue;
        else:
            for i in range(0,len(g1_pseudotime_locates)):
                gene_pseudotime_locates.append(g1_pseudotime_locates[i])
        
        if num%100 ==0:
            print('whole ',num,' genes have been done.')
    
    #func = partial(calculate_genes_pseudotime_location,rdata_df,cell_sets,markers_s_LGPS_rdata,adata_obsm_df);
    #gene_pseudotime_locate_dicts = pool.map(func,iter_genes);
    #pool.close();
    #pool.join();

    gene_pseudotime_locates_df = pd.DataFrame(gene_pseudotime_locates);
    
    end = datetime.datetime.now()
    #print(end);
    print('Running time : %s Seconds' %(end-start))
    print('Done!');
    return gene_pseudotime_locates_df;

def Fisher_test_for_each_gene(rdata_df: 'dataFrame of raw data',cell_sets:'all cell sets',num_allCells:'number of all cells',gene:'gene name'):
    """
    get ratio, p value and q value of Fisher test for given gene.
    
    Parameters
    ----------
    rdata_df : data frame format of raw annoData,
    cell_sets : a dict. key is the cell tyoe; value is a list of cell name.
    num_allCells: the number of all cells which is used as the background cells.
    gene: a given gene name
    
    return
    ----------
    gene_ratio : a dict of which key is gene and value is ratio 
    gene_pv : a dict of which key is gene and value is p value from fisher test
    gene_qv : a dict of which key is gene and claue is q value which is FDR for value by using BH method
    """
    gene = gene;
    gene_highly_cells = get_highly_cells_for_each_gene(rdata_df,gene);
    test_cells = gene_highly_cells[gene];
    num_test_cells = len(test_cells); 
    
    gene_pv = dict();
    gene_qv = dict();
    gene_ratio = dict();
    
    all_ratio = [];
    all_pv = [];
    all_qv = [];
    
    if num_test_cells == 0:
        gene_ratio[gene] = [];
        gene_pv[gene] = [];
        gene_qv[gene] = [];
        return (gene_ratio,gene_pv, gene_qv)
    else:
        for cs in cell_sets.keys():
            cell_set_cs = cell_sets[cs];
            num_cell_set = len(cell_set_cs);
            a = len(set(test_cells)&set(cell_set_cs));
            b = num_test_cells-a;
            c = num_cell_set-a;
            d = num_allCells-num_cell_set;
            ra = a/float(num_test_cells);
            pv = fisher.pvalue(a, b, c, d).right_tail;  # change
            #pv = fisher_exact([[a,b],[c,d]])[1];
            all_pv.append(pv);
            all_ratio.append(ra);
    
        all_qv = bh_qvalues(all_pv);
    
        gene_pv[gene] = all_pv;
        gene_qv[gene] = all_qv;
        gene_ratio[gene] = all_ratio;
        
        return (gene_ratio, gene_pv, gene_qv);

def bh_qvalues(pv):
    if pv == []:
        return [];
    m = len(pv);
    args,pv = zip(*sorted(enumerate(pv),key=operator.itemgetter(1)));
    if pv[0] < 0 or pv[-1] >1:
        raise ValueError("p-values must be between 0 and 1");
    qvalues = m*[0];
    mincoeff = pv[-1];
    qvalues[args[-1]] = mincoeff;
    for j in range(m-2,-1,-1):
        coeff = m*pv[j]/float(j+1);
        if coeff < mincoeff:
            mincoeff = coeff;
        qvalues[args[j]] = mincoeff;
    
    return qvalues
     
# get highly expressed cells for given gene
def get_highly_cells_for_each_gene(rdata_df: 'dataFrame of raw data', gene: 'gene name'):
    """
    get highly expressed cells for given gene.
    
    Parameters
    ---------
    rdata_df:
    gene : 
    """
    gene_list = rdata_df[gene]
    gene_filter = gene_list.loc[gene_list>0,]
    gene_sort_df_filter = gene_filter.sort_values(ascending=False)
    
    gene_highly_cells = dict();
    #gene_sort_df = rdata_df.sort_values(by=gene,ascending=False);
    #gene_sort_df_filter = gene_sort_df.loc[gene_sort_df[gene]>0,];
    # get the cells which highly express given gene
    if len(gene_sort_df_filter)<=10:
        #print('the '+gene+' has less than 10 highly expressed cells ');
        gene_highly_cells[gene] = []
        return gene_highly_cells
    
    else:
        matAA = curve_fitting(gene_sort_df_filter,gene);
        df_3_solution = calculate_derivation(matAA);
        inflection_point=0;
        for i in range(len(df_3_solution)):
            if str(df_3_solution[i])[-2] =='j':
                if str(df_3_solution[i])[-4:-1] =='+0j':
                    if df_3_solution[i] > 0:
                        inflection_point = math.ceil(np.real(df_3_solution[i]))
                        break;
                    else:
                        continue;
            else:
                if df_3_solution[i]>0:
                    inflection_point = math.ceil(np.real(df_3_solution[i]));
                    break;
        highly_cells = list(gene_sort_df_filter.index[0:inflection_point]);
    gene_highly_cells[gene] = highly_cells
    
    return gene_highly_cells


def calculate_genes_pseudotime_location(rdata_df: 'dataFrame of annoData',
                                        cell_sets,
                                        markers_s_LGPS_rdata,
                                        adata_obsm_df,
                                        gene: 'gene name'):
    gene = gene;
    gene_pseudotime_locate = list();
    
    gene_highly_cells = get_highly_cells_for_each_gene(rdata_df,gene)[gene];
    
    if gene_highly_cells ==[]:
        pass;
    else:
        if gene in list(markers_s_LGPS_rdata['gene_name']):
            gene_cell_types = list(markers_s_LGPS_rdata.loc[markers_s_LGPS_rdata['gene_name']==gene,'cell_type']);
            for CT in gene_cell_types:
                gene_pseudotime_locate_dict = dict();
                gene_cell_type_celllist = cell_sets[CT];
    
                gene_cell_type_highly_cell = list(set(gene_highly_cells) & set(gene_cell_type_celllist));
    
                gene_pseudotime_cells_df = adata_obsm_df.loc[adata_obsm_df.index.isin(gene_cell_type_highly_cell),:];
                gene_pseudotime = np.median(gene_pseudotime_cells_df['dpt_pseudotime']); 
                # calculating the x and y location for each gene
                gene_pseu_dist = abs(gene_pseudotime_cells_df['dpt_pseudotime']-gene_pseudotime);
                gene_min_pseu_dist_cell = gene_pseu_dist.idxmin();
        
                if len(gene_pseudotime_cells_df.columns) >= 4:
                    gene_x_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,0];
                    gene_y_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,1];
                    gene_z_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,2];
            
                    gene_pseudotime_locate_dict['gene_name'] = gene;
                    gene_pseudotime_locate_dict['x_location'] = gene_x_locate;
                    gene_pseudotime_locate_dict['y_location'] = gene_y_locate;
                    gene_pseudotime_locate_dict['z_location'] = gene_z_locate;
                    gene_pseudotime_locate_dict['dpt_pseudotime'] = gene_pseudotime;
        
                if len(gene_pseudotime_cells_df.columns) == 3:
                    gene_x_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,0];
                    gene_y_locate = gene_pseudotime_cells_df.loc[gene_min_pseu_dist_cell,1];
            
                    gene_pseudotime_locate_dict['gene_name'] = gene;
                    gene_pseudotime_locate_dict['x_location'] = gene_x_locate;
                    gene_pseudotime_locate_dict['y_location'] = gene_y_locate;
                    gene_pseudotime_locate_dict['dpt_pseudotime'] = gene_pseudotime;
                
                gene_pseudotime_locate.append(gene_pseudotime_locate_dict);
        else:
            pass;
    return gene_pseudotime_locate;

def projection(A,b):
    AA = A.T.dot(A)
    w = np.linalg.inv(AA).dot(A.T).dot(b)
    #print(w)
    return w

# curve fitting using ordinary least squares techniques
def curve_fitting(sort_data: 'sorted dataFrame of annoData',gene: 'gene name'):
    xa = np.arange(1,len(sort_data)+1,dtype=float)
    #xa=[x for x in range(1,len(sort_data[gene])+1)]
    ya=list(sort_data)
    #xa = np.array(xa)
    #xa.dtype="float64"
    y1 = np.array(ya)
    # power=9
    order=9;

    m = []

    for i in range(10):
        a = xa**(i)
        m.append(a)
    A = np.array(m).T
    b = y1.reshape(y1.shape[0],1)

    matAA = projection(A,b)
    matAA.shape = 10,

    return matAA
   
# calculating the first,second and third derivation of fitted curve
def calculate_derivation(matAA):
    df_3_diff = []
    for i in range(len(matAA)):
        df_3_diff.append(matAA[i]*i*(i-1)*(i-2))
    df_3_diff_reverse = list(reversed(df_3_diff[3:]))
    df_3_solution = np.roots(df_3_diff_reverse)
    df_3_solution.sort()
    return df_3_solution

def draw_fitted_curve(gene_sort_df_filter,gene):
    #draw the curve after fitting
    #print(matAA)
    xa=[x for x in range(1,len(gene_sort_df_filter[gene])+1)]
    ya=list(gene_sort_df_filter[gene])
    matAA = curve_fitting(gene_sort_df_filter,gene);
    xxa=[x for x in range(1,len(gene_sort_df_filter[gene])+1)]
    yya=[]
    for i in range(0,len(xxa)):
        yy=0.0
        for j in range(0,10):
            dy=1.0
            for k in range(0,j):
                dy*=xxa[i]
            dy*=matAA[j]
            yy+=dy
        yya.append(yy)
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(xa,ya,color='green',linestyle='',marker='o',markersize=4)
    ax.plot(xxa,yya,color='red',linestyle='-',marker='',linewidth=2)
    return fig
# Screening differentially expressed gene using this model