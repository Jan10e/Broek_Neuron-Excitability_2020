## Title: Calculation of coordinates of saddle-node (SN) bifurcation
#
# The saddle node is the situation where the V-nullcline and
# w-nullcline are tangential, irrespective of the epsilon value
############################################################

# Packages
using NLsolve
using DelimitedFiles

# Gating functions
winf1(v::Float64,v0::Float64,w0::Float64) = 10/(1 + exp(-5 * (v - (1/3) - v0))) + w0;

################ Ranges for convergence #################
const RANGE_LEFT = -1.1
const RANGE_RIGHT = -0.9

################### Get SN coorinates ##########################
# get v_sn coordinates when nullclines are tangential, i.e. det(J)=0
function FM_sn_coordinates(v0::Real, w0::Real, initClose::Real)

    v_sn = nothing
    w_sn = nothing
    I_app_sn = nothing

    # Get V-coordinate at SN
    detJSS(v) = (1-v^2)*(-1)-(-2*winf1(v,v0,w0))*((50*exp(5*v0 - 5*v + 5/3))/(exp(5*v0 - 5*v + 5/3) + 1)^2)
    detJSS(a::Array) = detJSS(a...) # "..." is the unpacking operator to unpack the array content

    sol = nothing
    loopCounter = 1
    while true
        sol = nlsolve(detJSS, [initClose]) # 0.88 seconds for the whole outer loop
        # println(sol)
        converged(sol) && break
        # converged(sol) returns true if successfully converged, in which case the break will be executed and exit the loop
        if(isnan(initClose))
            initClose = rand() * (RANGE_RIGHT-RANGE_LEFT) + RANGE_LEFT # random number between RANGE_LEFT and RANGE_RIGHT
        end
        if loopCounter >= 10
            break
        end
        loopCounter += 1
    end

    # Check whether algorithm converged
    if converged(sol) == true
            global v_sn = sol.zero[1]
    else
            global v_sn = NaN
    end

    # get n_sn and I_app_sn
    if isnan(v_sn) == false
        global w_sn = winf1(v_sn,v0,w0)
        global I_app_sn = -(v_sn - (1/3) * v_sn^3 - (10/(1 + exp(-5 *(v_sn - (1/3) - v0))) + w0)^2)
    else
        global w_sn = NaN
        global I_app_sn = NaN
    end

    return v_sn, w_sn, I_app_sn
end
