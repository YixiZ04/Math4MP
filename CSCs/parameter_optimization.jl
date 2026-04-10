
""" 
    This Script will perform a parameter optimization using Optuna, which used Tree-Structure Parzen Estimator (TPE) 
    as the optimization algorithm, an advanced Bayesian optimization method.
    
    4 parameters will be searched here if ran with no changing in the search space:
        - CSC growth rate
        - Progeny growth rate
        - Mature cell growth rate
        - Migration rate
    And 4 loss functions are being used here:
        - Difference to Iterations. The "objective" value has been set to 1800 (10 months with delta_t = 4h)
        - The Gompertzian growth curve fit (RMSE).
        - The population distribution for each cell type. The objective proportions are  [0.01, 0.04, 0.16, 0.158, 0.632]
        - The surface regularity ratio. The "objective" value here is 0.7.
    If not changing anything, 50 trials will be run and the search space will be:
        - CSC growth rate: [1, 500]
        - Progeny growth rate: [1, 150]
        - Mature cell growth rate: [100, 1000]
        - Migration rate: [80, 300]
    NOTE: If wanted to change search space, it should be changed at the objective functions found at the very end of this script.
    This outputs an .tsv file with no headers, but the values will be:
    CSC_rate    Progeny_rate    Mature_cell_rate    Mig_rate    RMSE    iter_diff   regularity_diff    proportions diff
""" 

using PythonCall
using Distributions     # Includes binomial and multinomial distributions
using Random            # Allows random sampling from previous probability distributions
using DelimitedFiles    # Enhances file I/O

include("constants.jl")
include("monitor.jl")
include("grid.jl")

optuna = pyimport("optuna")

# Gompertz growth model paramters

DELTAT = 4                                                                  # Time step length (hours)
K = 158.04                                                                  # Asymptotic volume (mL)
V0 = 1.0                                                                    # Initial volume (mL)
A = 0.01                                                                    # Gomperztian parameter (1/day)

# Objective parameters
MEAN_ITERATION = 1800                                                       # 10 * 30 * 24 / DELTAT. (10 months)
OBJECTIVE_SURFACE_REGULARITY = 0.7                                          # Other values could be set as well.
OBJECTIVE_POPULATION_PROPORTION = [0.01, 0.04, 0.16, 0.158, 0.632]          # Can alse be different.

# Optimization parameters

N_TRIALS = 1000                                                             # Change to 1000 to get a more robust search


# HELPER FUNCTIONS

function _get_last_non_zero_item(array::Array{Int64, 1})
    """ 
        This function will find and return the last non-zero item in an array.
        This is useful when trying to get the final population count for each cell type as there are many 0 values at the end of the array.
        To avoid any particular case, the entire array will be iterated.
    """
    non_zero_item = 0
    for item in array
        if item != 0
            non_zero_item = item
        else
            continue
        end
    end
    return non_zero_item

end


# GOMPERTZIAN CURVE FUNCTIONS
function _get_days_from_iteration(iterations::Vector{Int64}, deltat::Int64)
    """
        Assuming one day has 24 hours. From iteration number, this should return the number of days.    
    """
    return iterations * deltat / 24
end

function _get_gmb_gompertzian_curve(V0::Float64, K::Float64, a::Float64, t::Vector{Float64})
    """
        This function calculates the Gompertzian curve for Glioblastoma described by Stensjøen et al. (2015), which is given by the formula:
        V(t) = K * exp(log(V0/K) * exp(-a*t))
        Where V0 is the initial volume, K is the asymptotic volume, a is the Gompertzian parameter and t is the time in days.
    """
    return K .* exp.(log(V0/K) .* exp.(-a.*t))
end

# SURFACE REGULARITY FUNCTIONS
function _calculate_distance_voxels(i1::Int64, j1::Int64, k1::Int64, i2::Int64, j2::Int64, k2::Int64)
    """
        This function calculates the distance between two voxels with coordinates (i1,j1,k1) and (i2,j2,k2).
        This is useful to calculate the surface regularity, as we can compare the distance of each occupied voxel to the center with the radius of a sphere with the same volume as the tumor.
    """    
    return sqrt((i1 - i2)^2 + (j1 - j2)^2 + (k1 - k2)^2)
end

function _get_neighbors_count(occ::Array{CartesianIndex{3},1})
    """
        This functions will iterate through the array of occupied voxels and find the number of neighbors for each voxel.
        It outputs a 1D array with the neighbord count.
    """
    Neighbors_count = zeros(Int64, length(occ))
    voxel_count = 0
    for voxel1 in occ
        voxel_count += 1  
        for voxel2 in occ
            distance = _calculate_distance_voxels(voxel1[1], voxel1[2], voxel1[3], voxel2[1], voxel2[2], voxel2[3])
            if distance == 1
                Neighbors_count[voxel_count] += 1
            end
        end
    end
    return Neighbors_count
end

function _get_surface_by_occupied_voxels(Occ::Array{CartesianIndex{3},1})
    """
        This functions calculates the surface area from an array of occupied voxels: (x,y,z). 
        It is based on the fact that each voxel is a cube with area = 6 and we want to get the neighbours of each voxel.
        And the area of each voxel will be calculated as 6 - number of neighbours. So if a voxel has 6 neighbours, its area will be 0, and if it has no neighbours, its area will be 6.
    """
    surface_array = fill(6, length(Occ))
    neighbors_count = _get_neighbors_count(Occ)
    surface = sum(surface_array - neighbors_count)
    return surface
end

function _calculate_surface_regularity(Vol::Float64, Surface::Int64)
    """
        This function calculates a ratio Volume / Surface, assuming that the tumor is a perfect sphere. 
        The ratio = 6sqrt(pi)*Vol/Surface^(3/2).
        If ratio is close to 1, the tumor is regular, and if it is close to 0, the tumor is irregular.
    """
    return 6 * sqrt(pi) * Vol / Surface^(3/2)
end

# POPULATION FUNCTIONS

function _get_population_proportion(pops::Matrix{Float64})
    """
        This function calculates the proportion of each cell type in the total population.
        It outputs a 1D array with the proportion of each cell type.
    """ 
    subpopulation_count_array = [maximum(row) for row in eachrow(pops)]
    return subpopulation_count_array ./ sum(subpopulation_count_array)
end


# LOSS FUNCTIONS
function _iter_length_diff(mean_iteration::Int64, iterations::Int64)
    """ 
        Calculate the absolute difference between the mean number of iterations as the the survival time for Glioblastoma
        is between 8-12 months, so here the mean of 10 months has been considered, and assuming a deltat of 4 hours,
        the mean number of iterations would be 1800. 
        This function will return the absolute difference between this mean number of iterations and the number of iterations that the simulation has taken to reach the maximum volume, which is a value to minimize in the optimization process.
    """
    return abs(mean_iteration - iterations)
end

function _gompertz_rmse(iteration_array::Array{Int64,1}, simulation_array::Array{Float64,1})
    """
        This function calculates the RMSE between the simulated growth curve and the Gompertzian curve for Glioblastoma described by Stensjøen et al. (2015).
        This output is a single number and this will be the a valor to minimize in the optimization process. 
    """
    t = _get_days_from_iteration(iteration_array, DELTAT)
    gompertzian_curve = _get_gmb_gompertzian_curve(V0, K, A, t)
    rmse = sqrt(mean((simulation_array - gompertzian_curve).^2))
    return rmse
end

function _surface_regularity_diff(regularity_ratio::Float64, objective_regularity::Float64)
    """ 
        This function calculates the difference between the regularity ratio and a objective regularity ratio, which will be 0.7
    """
    return abs(regularity_ratio - objective_regularity)
end

function _proportion_diff(population_proportions::Array{Float64, 1}, objective_proportions::Array{Float64, 1})
    """ 
        This function gives back a the RMSE between the popolation proportions and the objective propotion, which will be:
        [0.01, 0.04, 0.16, 0.158, 0.632]
    """
    rmse = sqrt(mean((population_proportions - objective_proportions).^2))
    return rmse
end


# MAIN FUNCTION
function _main_process(g::Grid, c::Constants, m::Monitor)
    """
        This is be the main process function. 
        The 
    """
    iterations = Int64[]
    @time while m.Vol2[m.evalstep] < c.VolEnd

        grid_time_step!(g, c, m)    # Perform a single iteration

        increase_tstep(m)           # Move to the next time step

        # If current time step requires a system evaluation
        if m.t % c.NstepNevalRatio == 0
            # Update all monitored variables
            println("The iteration ", m.t ," is running...")
            println("The current volume: ", m.vol[m.evalstep], " mL")
            update_monitor_stats!(m, c)
            push!(iterations, m.t)
        end        
    end
    return (g.Occ, iterations, m.vol, m.pops)
end

# The objective function

function objective(trial)
    
    c = Constants()
    g = Grid(c)
    m = Monitor(c)

    # DEFINE THE SEARCH SPACE

    CSC_rate = pyconvert(Float64,trial.suggest_float("CSC_rate", 1, 500))
    Progeny_rate = pyconvert(Float64, trial.suggest_float("Progeny_rates", 1, 150))
    Mature_cell_rate = pyconvert(Float64, trial.suggest_float("Mature_cell_rate", 100, 1000))
    mig_rate = pyconvert(Float64, trial.suggest_float("mig_rate", 80, 300))
    c.Grate[1] = CSC_rate
    c.Grate[2] = c.Grate[3] = Progeny_rate
    c.Grate[4] = c.Grate[5] = Mature_cell_rate
    c.Migrate = fill(mig_rate, 5)

    Occ, iterations, volume, pops = _main_process(g, c, m)

    # Prepare for iteration diff loss function
    
    total_iteration = maximum(iterations)
    

    #Prepare for Gompertzian model fit
    iterations = push!(iterations, total_iteration + 20)
    no_zeros_volume_array = volume[1:length(iterations)]



    # Prepare for surface regularity loss function
    vol = maximum(volume)
    surface = _get_surface_by_occupied_voxels(Occ)
    regularity_ratio = _calculate_surface_regularity(vol, surface)

    # Prepare for proportion fit

    proportions = _get_population_proportion(pops)

    # Get losses
    rmse = _gompertz_rmse(iterations, no_zeros_volume_array)
    iter_loss = _iter_length_diff(MEAN_ITERATION,total_iteration)
    surface_regularity_loss = _surface_regularity_diff(regularity_ratio, OBJECTIVE_SURFACE_REGULARITY)
    proportion_loss = _proportion_diff(proportions, OBJECTIVE_POPULATION_PROPORTION)

    return rmse, iter_loss, surface_regularity_loss, proportion_loss
end

# Optimization process

study = optuna.create_study(
    directions = pylist(["minimize", "minimize","minimize","minimize"])
)
study.optimize(objective, n_trials=N_TRIALS)

# Get result files

trials = study.trials
results = []

for trial in trials

    CSC_rate = pyconvert(Float64, trial.params["CSC_rate"])
    Progeny_rate = pyconvert(Float64, trial.params["Progeny_rates"])    
    Mature_cell_rate = pyconvert(Float64, trial.params["Mature_cell_rate"])
    mig_rate = pyconvert(Float64, trial.params["mig_rate"])

    loss1 = pyconvert(Float64, trial.values[0])
    loss2 = pyconvert(Float64, trial.values[1])
    loss3 = pyconvert(Float64, trial.values[2])
    loss4 = pyconvert(Float64, trial.values[3])

    push!(results,
        [
            CSC_rate,
            Progeny_rate,
            Mature_cell_rate,
            mig_rate,
            loss1,
            loss2,
            loss3,
            loss4
        ]
    )

end

results_matrix = reduce(vcat, permutedims.(results))

writedlm(
    "optuna_results.tsv",
    results_matrix,
    '\t'
)