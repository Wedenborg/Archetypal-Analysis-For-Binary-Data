import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'



def plot_loss_arc(L,n_arc_list,i,color,model, dataset, savedir = None):


    fig, ax = plt.subplots(figsize=(15, 5),layout='constrained')

    ax.errorbar(n_arc_list,np.mean(L,axis=1),yerr=np.std(L,axis=1),c=color,label=f'{model} - {dataset}')

    ax.set_xticks(n_arc_list,fontsize=25)
    ax.set_xlabel('Number of Archetypes',fontsize=30)
    ax.set_ylabel('Loss',fontsize=30)

    plt.legend(fontsize=30)

    plt.yticks(fontsize=25)
    plt.xticks(fontsize=25)    
    #plt.title(f'{model} - {dataset}',fontsize=30)

    if savedir is not None:
        plt.savefig(savedir+f"/loss_layer_"+str(i)+".png")
        plt.close()
    else:
        plt.show()