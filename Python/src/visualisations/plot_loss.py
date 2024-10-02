import matplotlib.pyplot as plt
import numpy as np

import matplotlib

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'



def is_monotonic(x):
    dx = np.diff(x)
    return np.all(dx <= 0) or np.all(dx >= 0)

def plot_loss(title, loss, varexpl, offset=0):
    """
    Plots the loss curve and displays additional information.

    Args:
        title (str): Title of the plot.
        loss (ndarray): Array of loss values.
        varexpl (float): Variance explained by the model.

    Returns:
        None
    """

    # Create a new figure and axis object.
    fig, ax = plt.subplots()

    # Plot the loss curve.
    ax.plot(loss)

    # Set the x and y axis labels and title.
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Loss')
    ax.set_title(title)
    ax.xaxis.get_major_locator().set_params(integer=True)
    if offset != 0:
        ax.set_xticks(np.arange(len(loss)))
        ax.set_xticklabels(offset+np.arange(0, len(loss)))

    # Create a text string with additional information.
    textstr = '\n'.join((
        r'Is loss monotonic: %s' % (is_monotonic(loss), ),
        r'Variance explained: %.2f' % (varexpl, )))

    # Set the properties of the text box.
    props = dict(boxstyle='round', facecolor='steelblue', alpha=0.5)

    # Add the text box to the plot.
    ax.text(0.38, 0.96, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)

    # Show the plot.
    plt.show()