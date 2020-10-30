# Title: plots frequency-current curve for mirrored FHN
#
# Author: Jantine Broek
# Created: May 2019
######################################################

# working directory
cd("/Users/jantinebroek/Documents/03_projects/03_excitability/Broek_Neuron-Excitability_2020/src")

# Loads packages
using Plots
pyplot()
using DelimitedFiles
using LaTeXStrings
using Printf

# Function scripts
include("FM_sn_coordinates.jl")
include("FM_FI_curve.jl")

# Bifucation parameters
v0 = -0.26          # Fast negative conductance gain
w0 = -0.04          # Slow positive conductance gain
epsilon = 1e-2      # Time constant for w


############# Get crit values ###############
initClose = -1.0
(v_sn, n_sn, I_app_sn) = FM_sn_coordinates(v0, w0, initClose)


############# Input ###############
const dI = 1e-1
const I_app_start = 0.5
const I_app_stop = 10.
const I_app_range = collect(I_app_start:dI:I_app_stop)

str_input = "V$v0 w$w0 e$epsilon I_range$I_app_start:$dI:$I_app_stop"

############ Euler Solve Info ######################
const tend = (5/epsilon)*10
const dt = 1e-3
const t = collect(0:dt:tend)
const t500 = convert(Int64,50/dt) # Removes transient (first 500 ms)

# initial values
const v_init = -1.2
const w_init = w0

############ Get frequency ######################
@time (freq, v) = ficurve(v0, w0, I_app_range, epsilon, v_init, w_init, t)

############ Plot FI ######################
cd("../figures/")    #store in `figures` file

display(plot(I_app_range, freq,
    color = :black,
    linewidth = 4,
    xlabel = L"$I_{app}$",
    ylabel = L"freq ($Hz$)",
    ylims = (0,180),
    # title = "FI curve with V0=$V0, w0=$w0, ϵ=$teps \n
    #  f0 = $(@sprintf("%.4f", f0)), I0 = $I0",
    title = "FI curve with winf new2 V0=$v0, w0=$w0, ϵ=$epsilon,
    dI_ramp=$dI",
    titlefontsize = 11,
    legend = :none,
    # seriestype = :scatter,
    ticks = true,
    grid = false
    ))

savefig("FM_FI_" * str_input * ".eps")

############ Save Data ######################
cd("../data")   #store in `data` file
writedlm("FM_FI_" * str_input * ".dat", freq)
