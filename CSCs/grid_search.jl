"""
    This Julia Script is used to perform a grid search for 2 parameters:
        - CSC growth rate
        - Progeny growth rate
    The search space is defined at the very beginning of this Script, by default, the value for the search is:
        - CSC growth rate: [200.0, 500.0, 20]
        - Progeny growth rate: [1.0, 31.0, 5]
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
deltat_array = [2, 3, 4]
CSC_grate_array = range(100.0, 300.0, step=50)
Progeny_grate_array = range(1.0, 31.0, step=5)    


# Define the maximun months. The max_iteration number would depend to this value and the deltat value example, if deltat is 4, then max_iteration would be 15*30*24/4 = 2700.
max_months = 15

# Output file. Make the file. Overwrite the existing one if it exists.
output_file = joinpath(@__DIR__, "grid_search_results.txt")
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
        for deltat in deltat_array
            println("Running simulation for CSC growth rate = $CSC_grate, Progeny growth rate = $Progeny_grate, deltat = $deltat")
            num_iter = 0
            # Update constants with current parameter values
            c = Constants()
            c.Grate[1] = CSC_grate
            c.Grate[2] = Progeny_grate
            c.Grate[3] = Progeny_grate
            c.deltat = Int64(deltat)
            c.Nstep = floor(c.tspan / c.deltat)
            c.Neval = ceil(c.Nstep / 20) + 1
            c.NstepNevalRatio = round(c.Nstep / c.Neval)

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
                    print("Starting for the iteation ", num_iter, " with CSC growth rate = ", CSC_grate, " and Progeny growth rate = ", Progeny_grate, "\n")
                            update_monitor_stats!(m, c)

                    print_curr_stats(m)
                end
            end
            # Store results
            
            open(output_file, "a") do file
                println(file, "$CSC_grate $Progeny_grate $num_iter $deltat")
            end

            println("Completed: CSC growth rate = $CSC_grate, Progeny growth rate = $Progeny_grate, Iterations to reach volume threshold = $num_iter")
   
        end
    end
end

# Functions

function get_max_iteration(deltat, max_months)
    """ 
        This functions gets the maximun iteration number for a simulation based on its deltat and the maximum
        months. This avoids that excesively long simulations, which is not the most realistic cases.
    """
    return Int64(floor(max_months * 30 * 24 / deltat))
end