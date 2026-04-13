"""
    parameter_optimization_functions.py

    Contains functions to be used to analyse the parameter optimization results

"""

# IMPORT MODULES
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import gaussian_kde

# HELPER FUNCTIONS

def _kde_estimation(samples, df_column_name):
    """
        Plots the Gaussian Kernel estimation curve.
    """
    kde = gaussian_kde(samples)
    x = np.linspace(samples.min(), samples.max(), 1000)
    plt.plot(x,
             kde(x),
             label=f"KDE",
             color="black",
             linewidth=2
             )

# FUNCTIONS

def save_plot (path2plot):
    """
        This function saves the figure when calling it to the path indicated.
        The dpi will be 400
    """
    plt.savefig(path2plot, format='png', dpi=400)
    print(f"The figure has been saved as {path2plot}")


def make_plots(df, df_column_name, color, title):
    """
        This function makes the plots when the column is indicated.
        The plot is a histogram and a KDE estimation of the density distribution.
    """
    plt.figure(figsize=(7, 10))
    sample_values = df[df_column_name].values
    plt.hist(sample_values,
             bins="fd",
             edgecolor="black",
             density=True,
             color=color,
             alpha=0.5,
             label=f"Optimization results"
             )
    _kde_estimation(
        sample_values,
        df_column_name
    )

    plt.title(f"{title} posterior distribution estimation")
    plt.xlabel(df_column_name)
    plt.ylabel("Density")
    plt.xlim(sample_values.min(), sample_values.max())
    plt.legend(loc="best")


