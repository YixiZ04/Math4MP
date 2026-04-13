"""
    Performs analysis on parameter optimization results:
        - Density plot.
        - Kernel density estimation curve.
    And saves plots on a new directory.
"""

# IMPORT MODULES

import os
import sys

import pandas as pd

from CSCs.parameter_optimization_analysis.parameters_analysis_functions import *

# DEFINE PARAMETERS
ROOT_DIR = os.getcwd()
PATH2INPUT = os.path.join(ROOT_DIR, "optimization_results.tsv")
PATH2OUTPUT = os.path.join(ROOT_DIR, "optimization_results_plots/")

# The main process

if __name__ == "__main__":
    print(f"Trying to fetch data from {PATH2INPUT}...")
    try:
        df = pd.read_csv(PATH2INPUT, sep="\t", dtype=np.float64, encoding="utf-8")
        print ("THE input data correctly fetched !")
    except FileNotFoundError:
        raise FileNotFoundError()

    print (f"Making the output directory{PATH2OUTPUT}...")
    os.makedirs(PATH2OUTPUT, exist_ok=True)

    print("Making plots for CSC rate...")
    make_plots(df,
               "CSC_rate",
               "green",
               "CSC growth rate"
               )
    save_plot(f"{PATH2OUTPUT}CSC_rate.png")
    plt.show()

    print("Making plots for Progeny rate...")
    make_plots(df,
               "Progeny_rates",
               "blue",
               "Progeny growth rate"
               )
    save_plot(f"{PATH2OUTPUT}Progeny_rates.png")
    plt.show()

    print("Making plots for Mature cell rates...")
    make_plots(df,
               "Mature_cell_rates",
               "purple",
               "Completly differentiated cell growth rate"
               )
    save_plot(f"{PATH2OUTPUT}Mature_cell_rates.png")
    plt.show()

    print("Making plots for migration_rate...")
    make_plots(df,
               "Migration_rate",
               "yellow",
               "Migration rate"
               )
    save_plot(f"{PATH2OUTPUT}Migration_rate.png")
    plt.show()

    print(f"All the plots have been successfully made! They are saved in {PATH2OUTPUT}.")
    print("Exiting the program with code 0...")
    sys.exit(0)



