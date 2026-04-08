"""
    This Julia Script is used to perform a grid search for 2 parameters:
        - CSC growth rate
        - Progeny growth rate
        - deltat (time step length)
    The search space is defined at the very beginning of this Script, by default, the value for the search is:
        - CSC growth rate: [100.0, 300.0, step=50]
        - Progeny growth rate: [1.0, 31.0, step=5]
        - deltat: [2, 3, 4]
    A limit has been set, as 15 months use to be the maximun survival time for patients diagnosed with Glioblstoma Multiforme.
    The maximun number of iterations will be calculated accorrding to the deltat values, assuming that:
        - 1 month = 30 days
        - 1 day = 24 hours
    So, if a simulations iterates for more than the maximun number of iterations, it will be stopped to avoid excessively long simulations,
    which are not the most realistic cases.rm 
    Usage: run on terminal julia grid_search.jl and the search will start. And the result file is "grid_search_results.txt" in the same folder as this script.
    The result analysis will be done in a separated script.

"""


using Distributions     # Includes binomial and multinomial distributions
using Random            # Allows random sampling from previous probability distributions
using DelimitedFiles    # Enhances file I/O

# Include files

include("constants.jl")
include("tools.jl")   
include("grid.jl")

# Define grid search ranges
CSC_grate_array = range(100.0, 200.0, step=10)
Progeny_grate_array = range(1.0, 11.0, step=1)    


# Define the maximun months. The max_iteration number would depend to this value and the deltat value example, if deltat is 4, then max_iteration would be 15*30*24/4 = 2700.
max_months = 15
deltat = 4

# Output file. Make the file. Overwrite the existing one if it exists.
output_file = joinpath(@__DIR__, "grid_search_results_2.txt")
if isfile(output_file)
    rm(output_file)
    touch(output_file)
else
    touch(output_file)
end

# Results storage   
println("Starting grid search for CSC growth rate and Progeny growth rate...")
for CSC_grate in CSC_grate_array
    for Progeny_grate in Progeny_grate_array
        println("Running simulation for CSC growth rate = $CSC_grate, Progeny growth rate = $Progeny_grate")
        num_iter = 0
        # Update constants with current parameter values
        c = Constants()
        c.Grate[1] = CSC_grate
        c.Grate[2] = Progeny_grate
        c.Grate[3] = Progeny_grate
    
        # Get the maximum iteration number based on deltat and max_months
        max_iteration = get_max_iteration(deltat, max_months)

        # Run simulation
        g = Grid(c)
        m = Monitor(c)

        while m.Vol2[m.evalstep] < c.VolEnd && num_iter < max_iteration
            grid_time_step!(g, c, m)
            increase_tstep(m)
            num_iter += 1
            if m.t % c.NstepNevalRatio == 0
                update_monitor_stats!(m, c)
                print_curr_stats(m)
            end
        end
        # Store results
        
        open(output_file, "a") do file
            println(file, "$CSC_grate $Progeny_grate $num_iter")
        end

        println("Completed: CSC growth rate = $CSC_grate, Progeny growth rate = $Progeny_grate, Iterations to reach volume threshold = $num_iter")

        
    end
end

# Functions

function get_max_iteration(deltat, max_months)
    """ 
        This functions gets the maximun iteration number for a simulation based on its deltat and the maximum
        months. This avoids that excesively long simulations, which is not the most realistic cases.
    """
    return Int64(max_months * 30 * 24 / deltat)
end
