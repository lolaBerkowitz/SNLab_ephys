import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
# Read the existing CSV file
os.chdir(os.path.dirname(__file__))
data = pd.read_csv('Combined.csv')
pivot = np.round(pd.pivot_table(data, 
                                index='group', 
                                aggfunc=np.mean),2)
numeric_columns = pivot.select_dtypes(include=[float]).columns
for col in numeric_columns:
    pivot[col] = (pivot[col]*100).round()
#If you want barplot, uncomment line below
#pivot.plot.barh(figsize=(10,7),title='Processing Progress')
plt.figure(figsize=(15,12))
heatmap=sns.heatmap(pivot,annot=True,cmap='YlGnBu_r',fmt=".2f",cbar_kws={'label':'Percent Completion'})
heatmap.set_xticklabels(heatmap.get_xticklabels(),rotation=45,ha='right')
heatmap.set_yticklabels(heatmap.get_yticklabels(),rotation=0,ha='right')
heatmap.figure.axes[1].yaxis.label.set_size(17)
plt.subplots_adjust(bottom=0.2)
plt.title('Processing Progress',y=1.05,fontdict={'size':30})
heatmap.set_ylabel('Group',fontdict={'size':17})
#TODO Change /Progress to \Progress for Windows
plt.savefig(os.getcwd()+'\Progress.png')
plt.show()