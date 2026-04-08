
""" 
    This Script will perform a parameter optimization using Optuna, which used Tree-Structure Parzen Estimator (TPE) 
    as the optimization algorithm, an advanced Bayesian optimization method.
    
    Many parameters should be tuned here, as the growth rates for each cell type. However, the main difference to the
    grid_search.jl is that here more results are taking into account, no only the time to reach the maximum volume,
    but also:
    - The Gompertzian growth curve. RMSE will be used as val loss.
    - The popuulation distribution for each cell type.
    - The surface regularity, which is quantified using the ratio between the tumor's volume and the volume of a sphere
    with the same surface area as the tumor.
""" 

using Optuna
using Distributions     # Includes binomial and multinomial distributions
using Random            # Allows random sampling from previous probability distributions
using DelimitedFiles    # Enhances file I/O

include("constants.jl")
include("grid.jl")
include("monitor.jl")

# Parameters considered for this process.

DELTAT = 4                         # Time step length (hours)
MEAN_ITERATION = 1800              # 10 * 30 * 24 / DELTAT. (10 months)
K = 158.04                         # Asymptotic volume (mL)
V0 = 0.001                         # Initial volume (mL)
A = 0.007545                       # Gomperztian parameter (1/day)

# Helper functions are built here to help the main process

function _get_days_from_iteration(iterations::Array{Int64,1}, deltat::Int64)
    """
        Assuming one day has 24 hours. From iteration number, this should return the number of days.    
    """
    return iterations * deltat / 24
end

function _get_gmb_gompertzian_curve(V0::Float64, K::Float64, a::Float64, t::Array{Float64,1})
    """
        This functions returns the Gompertzian curve for Glioblastoma described by Stensjøen et al. (2015).
        The formula is the following:
            V(t) = K*exp(log(V0/K)*exp(-a*t))
        V0: Initial volume (mL)
        K: asymptotic volume (mL)
        a: growth rate (1/day)
        t: time (days)
    """
    return K*exp(log(V0/K)*exp(-a*t))
end

# Define functions for losses.
function _calculate_difference_iter_length(mean_iteration::Int64=MEAN_ITERATION, iterations::Int64)
    """ 
        Calculate the absolute difference between the mean number of iterations as the the survival time for Glioblastoma
        is between 8-12 months, so here the mean of 10 months has been considered, and assuming a deltat of 4 hours,
        the mean number of iterations would be 1800. 
        This function will return the absolute difference between this mean number of iterations and the number of iterations that the simulation has taken to reach the maximum volume, which is a value to minimize in the optimization process.
    """
    return abs(mean_iteration - iterations)
end

function _calculate_rmse(iteration_array::Array{Float64,1}, simulation_array::Array{Float64,1})
    """
        This function calculates the RMSE between the simulated growth curve and the Gompertzian curve for Glioblastoma described by Stensjøen et al. (2015).
        This output is a single number and this will be the a valor to minimize in the optimization process. 
    """
    t = _get_days_from_iteration(iteration_array, DELTAT)
    gompertzian_curve = _get_gmb_gompertzian_curve(V0, K, A, t)
    rmse = sqrt(mean((simulation_array - gompertzian_curve).^2))
    return rmse
end


# The objective function

function objective(trial::Trial)

end

