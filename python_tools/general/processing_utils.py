import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from neuro_py.plotting.figure_helpers import set_size
import seaborn as sns

def plot_data_progress(summary_table):

    fig = plt.figure(figsize=set_size('thesis', fraction=1.5, subplots=(1, 1)))
    fig.subplots_adjust(bottom=0.2)

    heatmap = sns.heatmap(
        summary_table,
        annot=True,
        cmap="YlGnBu_r",
        fmt=".2f",
        cbar_kws={"label": "Percent Completion"},
    )

    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, ha="right")
    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, ha="right")
    heatmap.figure.axes[1].yaxis.label.set_size(17)
    plt.title("Processing Progress", y=1.05, fontdict={"size": 30})
    heatmap.set_ylabel("Group", fontdict={"size": 17})
    # TODO Change /Progress to \Progress for Windows
    plt.savefig(os.getcwd() + os.sep +  "Progress.png")
    plt.show()

def load_data_summary(path_to_csv):

    # import data for subject ID & basepath
    data = pd.read_csv(os.path.join(path_to_csv,"Combined.csv"))

    # get mean, as that will tell us propotion processed
    summary_table = np.round(pd.pivot_table(data, index="group", aggfunc=np.mean), 2)
    numeric_columns = summary_table.select_dtypes(include=[float]).columns

    # multiply by 100 to get percentage
    for col in numeric_columns:
        summary_table[col] = (summary_table[col] * 100).round()

    
    return summary_table