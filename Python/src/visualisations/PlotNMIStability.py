import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def plot_nmi_stability(NMI,n_arc_list,model,dataset, colors, savedir = None):
    
    
    df = pd.DataFrame(NMI.T)
    df['Method'] = f'Model {model}'


    df = df.melt(id_vars='Method', var_name='Archetypes', value_name='NMI')
    fig, ax = plt.subplots(1,1,figsize = (15,5), layout='constrained')

    ax = sns.boxplot(x='Archetypes', y="NMI", hue="Method", showmeans=True, data=df,palette=colors,meanprops={"marker": "s", "markerfacecolor": "white", "markeredgecolor": "black"})
    ax.xaxis.grid(True, which='major')
    [ax.axvline(x+.5,color='k') for x in ax.get_xticks()]
    plt.xticks(np.arange(len(n_arc_list)),n_arc_list,fontsize=25)
    plt.yticks(fontsize=25)    
    ax.set_xlabel('Number of archetypes', fontsize=30)
    ax.set_ylabel('NMI', fontsize=30)

    ax.set_ylim([0,1.05])
    plt.legend(fontsize=30,loc = 'lower right')
    plt.title(f'{model} - {dataset}',fontsize=30)
    
    if savedir is not None:
        plt.savefig(savedir+f"/NMI_model_"+str(model)+"dataset_"+str(dataset)+".png")
        plt.close()
    else:
        plt.show()

