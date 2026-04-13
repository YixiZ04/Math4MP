"""
    Paramters have been changed for this model:
        CSC_rate = [140, 160]
        Progeny_rate = [1, 10]
        Mature_cell_rate = [900, 1000]
        Mig_rate = [175, 210]    
        Death_rate = [80, 400]   (Only applied for the mature cells)
    This also implies that the "param_dist.txt" has also been changed adapting to these changes.
    
    Datatypes have also been changed:
        Grate::Array(Float64, 1)
        Drate::Array(Float64, 1)
        Migrate::Array(Float64, 1)
    Eliminated Mutrate, as it is not relevant and also will add an extra complexity to this model based on the CSC hypothesis.
"""



mutable struct Constants

    ################################################################################
    # SPATIO-TEMPORAL DOMAIN PARAMETERS
    ################################################################################

    TimeStart::Float64          # Use it to keep track of execution time
    deltat::Int64               # Time step length (hours)
    tspan::Int64                # Max number of simulation time allowed (hours)
    Nstep::Int64                # Max number of steps allowed
    N::Int64                    # Grid size (number of voxels per dimension: N x N x N)
    Neval::Int64                # Number of time steps between system evaluations
    NstepNevalRatio::Float64    # Number of system evaluations during whole simulation
    VolEnd::Float64             # Maximum tumor volume allowed. Once a tumor reaches this volume, simulation is stopped

    ################################################################################
    # STARTING CONDITIONS
    ################################################################################

    alt::Int64          # Number of alterations (number of clonal populations will be 2^alt)
    P0::Float64         # Initial cell number
    K::Int64            # Voxel's carrying capacity
    threshold::Float64  # Number of cells to be exceeded in order to consider a voxel active


    ################################################################################
    # CHARACTERISTIC TIMES
    ################################################################################

    fdata::Array{Float64,2}         # Input file containinig mean and standard deviation of cell processes' characteristic times
    CSC_mean::Float64               # Mean cell division time
    CSC_sd::Float64                 # Standard deviation of cell division time
    Progeny_mean::Float64
    Progeny_sd::Float64
    Matcell_mean::Float64
    Matcell_sd::Float64
    Drate_mean::Float64             # Mean cell death time
    Drate_sd::Float64               # Standard deviation of cell death time
    Migrate_mean::Float64           # Mean cell migration time
    Migrate_sd::Float64             # Standard deviation of cell migration time
    Grate::Array{Float64,1}         # Basal cell division time
    Drate::Array{Float64,1}         # Basal cell death time
    Migrate::Array{Float64,1}       # Basal cell migration time


    ################################################################################
    # AUXILIARY KERNEL FOR DISTRIBUTING MIGRATING CELLS IN THE NEIGHBOURHOOD
    ################################################################################

    c_old::Int64                # Iterator
    wcube::Array{Float64, 1}    # 3x3x3 kernel containing normalized migration probabilities to neighbour voxels, depending on distance from central voxel

    Pasim::Array{Float64, 1}
    Pchoice::Array{Float64, 2}

    ################################################################################
    # INITIALIZATION
    ################################################################################

    function Constants()

        TimeStart = time()

        deltat = 4
        tspan = 1e5
        Nstep = floor(tspan / deltat)
        N = 80
        Neval = ceil(Nstep / 20) + 1
        NstepNevalRatio = round(Nstep / Neval)
        VolEnd = 1.2e5;
        alt = 2
        P0 = 1e1
        K = 2e5
        threshold = 0.2 * K

        # Retrieve parameters from input file
        fdata = readdlm(joinpath(@__DIR__,"Param_dist.txt"))
        CSC_mean = fdata[1,1];                                  
        CSC_sd = fdata[1,2];
        Progeny_mean = fdata[2, 1];
        Progeny_sd = fdata[2, 2];
        Matcell_mean = fdata[3, 1];
        Matcell_sd = fdata[3, 2];
        Migrate_mean = fdata[4,1];
        Migrate_sd = fdata[4,2];
        Drate_mean = fdata[5,1];
        Drate_sd = fdata[5,2];
     
 
        # Random sample characteristic times from uniform distributions based in Param_dist.txt data
        Grate = [
                 rand(Uniform(CSC_mean-CSC_sd, CSC_mean+CSC_sd)),
                 rand(Uniform(Progeny_mean-Progeny_sd, Progeny_sd)),
                 rand(Uniform(Progeny_mean-Progeny_sd, Progeny_sd)),
                 rand(Uniform(Matcell_mean-Matcell_sd, Matcell_mean+Matcell_sd)),
                 rand(Uniform(Matcell_mean-Matcell_sd, Matcell_mean+Matcell_sd))
                ] 
        Migrate = rand(Uniform(Migrate_mean-Migrate_sd, Migrate_mean+Migrate_sd), 5)
        Drate = [0.0, 
                 0.0, 
                 0.0, 
                 rand(Uniform(Drate_mean-Drate_sd, Drate_mean+Drate_sd)), 
                 rand(Uniform(Drate_mean-Drate_sd, Drate_mean+Drate_sd))
                ]

        # Create weights for surrounding voxels (Moore neighbourhood)
        c_old = 0
        wcube = zeros(26)
        sumcube = 0

        # Create auxiliary kernel for distribution of migrating cells
        for i in [-1, 0, 1]
            for j in [-1, 0, 1]
                for k in [-1, 0, 1]
                    if abs(i) + abs(j) + abs(k) != 0
                        c_old = c_old + 1
                        wcube[c_old] = 1 / sqrt(abs(i) + abs(j) + abs(k))
                        sumcube = sumcube + wcube[c_old]
                    end
                end
            end
        end

        c_old = 0
        for i in [-1, 0, 1]
            for j in [-1, 0, 1]
                for k in [-1, 0, 1]
                    if abs(i) + abs(j) + abs(k) != 0
                        c_old = c_old + 1
                        wcube[c_old] = wcube[c_old] / sumcube
                    end
                end
            end
        end

        Pasim = [0.95, 0.95, 0.95, 0, 0]
        # Pchoice = [0:0:0:0:0 0.3:0:0:0:0 0.7:0:0:0:0 0:1:0:0:0 0:0:1:0:0]
        Pchoice = zeros(5, 5)
        Pchoice[1, 2] = 0.3
        Pchoice[1, 3] = 0.7
        Pchoice[2, 4] = 1
        Pchoice[3, 5] = 1

        new(TimeStart, deltat, tspan, Nstep, N, Neval, NstepNevalRatio, VolEnd, alt, P0, K, threshold, fdata,
        CSC_mean, CSC_sd, 
        Progeny_mean, Progeny_sd,
        Matcell_mean, Matcell_sd,
        Migrate_mean, Migrate_sd, 
        Drate_mean, Drate_sd, 
        Grate, Drate, Migrate, c_old, wcube, Pasim, Pchoice)

    end
end

