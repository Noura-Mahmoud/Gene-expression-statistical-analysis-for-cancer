import pandas as pd # to convert our data to dataframes
import matplotlib.pyplot as plt # to plot the correlation coefficient

#opening our data files

healthy = pd.read_csv('Desktop/Statistics/lusc-rsem-fpkm-tcga_paired.txt', sep='\t') # reading The file of healthy genes expressions data , Files are tab-separated.
cancer = pd.read_csv('Desktop/Statistics/lusc-rsem-fpkm-tcga-t_paired.txt', sep='\t') # reading The file of healthy genes expressions data , Files are tab-separated.

pd.options.display.max_columns = None

# Correlation Coefficient and Plotting

# filtration of healthy dataframe by dropping rows that contains more than 25 zeros in columns(genes expressions) 
fhealthy = healthy[healthy.astype(bool).sum(1)>27] #this func counts the 0's in each row and if it's >= 25 it returns false in the row else it returns true and then 
# we drop the rows that returned false from the original dataframe and filtration is done of data that gives results that doesn't make sense  
fcancer = cancer[healthy.astype(bool).sum(1)>27]

fhealthy = fhealthy[cancer.astype(bool).sum(1)>27]
fcancer = fcancer[cancer.astype(bool).sum(1)>27]

# computing correlation Coefficient for each row and its opposite in the 2 dataframes  (healthy , cancer) and storing them in a list
CCo = []

from scipy.stats import pearsonr # func to determine pearson correlation coefficient 
for i in range (len(fhealthy)):  # this loop iterates over all the rows in both dataframes and counts r for every gene in the 2 dataframes 
    G_h = fhealthy.iloc[i, 2:]
    G_c = fcancer.iloc[i, 2:]
    r, _ = pearsonr(G_h, G_c)
    CCo.append(r)

fhealthy.insert(2, "CCo", CCo , True)  # inserting column of CCo list to the filtered dataframe
ranked_healthy = fhealthy.sort_values('CCo')  # ranking this dataframe according to the values of CCo


Max_CCo = max(CCo) #getting the highest positive r 
Min_CCo = min(CCo) #getting the lowest negative r 
max_index = CCo.index(Max_CCo)  #getting the index of max r
min_index = CCo.index(Min_CCo)  #getting the index of min r

print(max_index)
print(min_index)

Max_CCo_gene = fhealthy.iloc[ max_index , : ] #getting the gene that has highest positive r 
Min_CCo_gene = fhealthy.iloc[ min_index , : ] #getting the gene that has lowest negative r 

plt.scatter(fhealthy.iloc[max_index, 2:] , fcancer.iloc[10790, 2:])  # Plotting the gene expressions of healthy and cancer gene that has Max r 
plt.title('Max Correlation Coefficient Plot') #giving title and labels specific names in the plot
plt.xlabel('Healthy')
plt.ylabel('Cancer')
plt.show()

plt.scatter(fhealthy.iloc[min_index, 2:] , fcancer.iloc[12926, 2:]) # Plotting the gene expressions of healthy and cancer gene that has Min r 
plt.title('Min Correlation Coefficient Plot') #giving title and labels specific names in the plot
plt.xlabel('Healthy')
plt.ylabel('Cancer')
plt.show()

