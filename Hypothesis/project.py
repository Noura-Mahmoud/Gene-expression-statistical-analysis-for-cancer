import pandas as pd
from scipy.stats import ttest_ind, ttest_rel
from statsmodels.stats.multitest import multipletests


healthy = pd.read_csv('data/lusc-rsem-fpkm-tcga_paired.txt', sep='\t')
cancer = pd.read_csv('data/lusc-rsem-fpkm-tcga-t_paired.txt', sep='\t')
pd.options.display.max_columns = None

healthy1 = healthy[healthy.astype(bool).sum(1)>27]
cancer1 = cancer[healthy.astype(bool).sum(1)>27]

healthy2 = healthy1[cancer1.astype(bool).sum(1)>27]
cancer2 = cancer1[cancer1.astype(bool).sum(1)>27]

#new length of arrays is 17275

#for independent
accp_i=[] #accepted probability list
accpgenes_i=[] #rejected probability genes list
rep_i=[] #accepted probability list
repgenes_i=[]  #rejected probability genes list

#for paired
accp_r=[] #accepted probability list
accpgenes_r=[] #accepted probability genes list
rep_r=[] #rejected probability list
repgenes_r=[]  #rejected probability genes list

p_values_r=[] # propability list for paired
p_values_i=[]  # propability list for independent 

#Assume confidence level=95%, alpha=0.05

for m in range(len(healthy2)):
    G_h = healthy2.iloc[m, 2:] # piking gene no. m from healthy set
    G_c = cancer2.iloc[m, 2:] # piking gene no. m from cancer set

    p_val_r = ttest_rel(G_h, G_c).pvalue #paired
    p_values_r.append(p_val_r)  # propability list while paired

    if p_val_r >= 0.05:  #accept H0 as p > 0.05   # for paired
        accp_r.append(p_val_r)
        accpgenes_r.append(healthy2.iloc[m, 0]) 

    else:   #rejection region --> reject H0  
        rep_r.append(p_val_r)
        repgenes_r.append(healthy2.iloc[m, 0])

    p_val_i = ttest_ind(G_h, G_c).pvalue #independent
    p_values_i.append(p_val_i)  # propability list while independent

    if p_val_i >= 0.05:  #accept H0 as p > 0.05  #for independent 
        accp_i.append(p_val_i)
        accpgenes_i.append(healthy2.iloc[m, 0])

    else:   #rejection region --> reject H0  
        rep_i.append(p_val_i)
        repgenes_i.append(healthy2.iloc[m, 0])

#2a

print ('There are '+str (len(accp_i))+' accepted genes as independent')
# print (accpgenes_i)
print ('There are '+str (len(rep_i))+' rejected genes as independent')
# print (repgenes_i)
print ('There are '+str (len(accp_r))+' accepted genes as paired') 
# print (accpgenes_r)
print ('There are '+str (len(rep_r))+' rejected genes as paired')
# print (repgenes_r) 

##correction

#2b

#for independent
corrected_pi_values = multipletests(p_values_i, alpha=0.05, method='fdr_bh')[1]

#for paired
corrected_pr_values = multipletests(p_values_r, alpha=0.05, method='fdr_bh')[1]

# print('Independent corrected values : ', corrected_pi_values)
# print('Paired corrected values : ', corrected_pr_values)

significance_genes_i = pd.DataFrame({'p-values-independent':p_values_i, 'p-values_i_fdr':corrected_pi_values})
significance_genes_r = pd.DataFrame({'p-values-paired':p_values_r, 'p-values_r_fdr':corrected_pr_values})

#2c

#after applying correction method, we need to know the genes which has p<0.05 after correction
corrected_rejection_i = [] # here, we make a new list of corrected independent p values whose p=0.05, that will reject H0
corrected_rejection_r = [] # here, we make a new list of corrected paired p values whose p<0.05, that will reject H0

genes_corrected_rejection_i = [] # here, we make a new 'gene' list of corrected independent values whose p<0.05, that will reject H0
genes_corrected_rejection_r = [] # here, we make a new 'gene' list of corrected paired values whose p<0.05, that will reject H0


for q in range(len(healthy2)):

    if corrected_pr_values[q] < 0.05:
        corrected_rejection_r.append(corrected_pr_values[q])
        genes_corrected_rejection_r.append(healthy2.iloc[q, 0])

    if corrected_pi_values[q] < 0.05:
        corrected_rejection_i.append(corrected_pi_values[q])
        genes_corrected_rejection_i.append(healthy2.iloc[q, 0])

#priniting lists of independent and paired genes in rejection region after correction
# print('Independent corrected rejection genes : ', genes_corrected_rejection_i)
# print('Paired corrected rejection genes : ', genes_corrected_rejection_r)

significance_genes_i['significance:p_vlaue_i'] = significance_genes_i['p-values-independent'].apply(lambda x: x < 0.05)
significance_genes_i['significance:p_vlaue_i_fdr'] = significance_genes_i['p-values_i_fdr'].apply(lambda x: x < 0.05)
# Get significant genes after fdr correction
diffrentially_genes_i = significance_genes_i[significance_genes_i['significance:p_vlaue_i_fdr']== True]
# diffrentially_genes_i

significance_genes_r['significance:p_vlaue_r'] = significance_genes_r['p-values-paired'].apply(lambda x: x < 0.05)
significance_genes_r['significance:p_vlaue_r_fdr'] = significance_genes_r['p-values_r_fdr'].apply(lambda x: x < 0.05)
# Get significant genes after fdr correction
diffrentially_genes_r = significance_genes_r[significance_genes_r['significance:p_vlaue_r_fdr']== True]
# diffrentially_genes_r

#2d

genes_set_r = set(genes_corrected_rejection_r) # list conversion to set to apply builtin fns of sets
genes_set_i = set(genes_corrected_rejection_i)

intersect = genes_set_r.intersection(genes_set_i) # intersected genes
diff_ri = genes_set_r.difference(genes_set_i)  # r diff from i
diff_ir = genes_set_i.difference(genes_set_r)  # i diff from r

print('independent intersected with paired = ', len(intersect))
print('paired - independent = ', len(diff_ri))
print('independent - paired = ', len(diff_ir))