## Title: simulation of mirrored FHN model using
# the step-down protocol.
#
# Author: Jantine Broek
# Date: May 2019
##################################################

# working directory
cd("/Users/jantinebroek/Documents/03_projects/03_excitability/Broek_Neuron-Excitability_2020/src")

# Loads packages
using Plots     # Select the Plots package for figures
pyplot()        # Plots package will use Pyplot (matplotlib needs to be installed).
using LinearAlgebra
using DelimitedFiles
using LaTeXStrings
using Printf

# Function scripts
include("FM_FI_curve_step_down.jl")

# Functions
winf1(v::Real) = 10/(1 + exp(-5 * (v - (1/3))))
winf3(v::Float64,v0::Float64,w0::Float64) = 10/(1 + exp(-5 * (v - (1/3) - v0))) + w0

############ Model parameters ######################
v0 = 0.64           # Fast negative conductance gain
w0 = -0.05          # Slow positive conductance gain
epsilon = 1e-2      # Time constant for w

############ Euler Solve Info ######################
# Simulation parameters
const T = round(5/epsilon)*10
const dt = 1e-3
const Tdt = convert(Int64,T/dt)
const t = range(dt,stop=T,length=Tdt)

# initial values
const v_init = 0. #0 for step-down, i.e. bistability; -1.2, when used for step up
const w_init = w0


############# Input ###############
const dI_step = 1e-1
const I_step_start = 0.5
const I_step_stop = 10.
const I_step_range = collect(I_step_start:dI_step:I_step_stop)

# steps down is at 1/5 of T, T500 is at 3/5 of T
const T_step_start = T/5
const T_step_stop = T
const I_app = 10.
const Ttrans = convert(Int64,((3*T)/5)/dt) # removes non-step period and transient

str_input = "V$v0 w$w0 e$epsilon I_range$I_step_start:$dI_step:$I_step_stop"

############ Get frequency ######################
@time (freq, v) = ficurve_step_down(v0, w0, epsilon, I_app, I_step_range)


############ Plot FI & Save ######################
cd("../figures/")    #store in `figures` file

display(plot(I_app .- I_step_range, freq,
    color = :black,
    linewidth = 2,
    xlabel = L"$I_{app}$",
    ylabel = L"freq ($Hz$)",
    # ylims = (0,12),
    # title = "FI curve with V0=$V0, w0=$w0, ϵ=$teps \n
    #  f0 = $(@sprintf("%.4f", f0)), I0 = $I0",
    title = "FI curve V0=$v0, w0=$w0, ϵ=$epsilon,
    dI_step=$dI_step",
    titlefontsize = 11,
    legend = :none,
    seriestype = :scatter,
    ticks = true,
    grid = false
    ))

savefig("FM_FI_step_down_" * str_input * ".eps")

############ Save Data ######################
cd("../data")   #store in `data` file
writedlm("FM_FI_step_down_" * str_input * ".dat", freq)
