## Title: Get multiple fixed point in mFHN model
#
# Author: Jantine Broek
# Date: Nov 2019
##################################################

# Packages
using NLsolve

function FM_fixed_points_multi(steady_state)

    # get all fixed points
    x_init = collect(-3.0:1e-2:3.0)
    fps = zeros(length(x_init))
    for i = 1:length(x_init)
        fps_temp = nlsolve(steady_state, [x_init[i]], autodiff = :forward, xtol = 1e-14, ftol = 1e-14, iterations = 1_00)
        global fps[i] = fps_temp.zero[1] #zero contains solution that has converged
    end

    # extract similar fixed points
    fps_rounded = Array{Real, 1}(undef, length(fps))
    for j = 1:length(fps)
        global fps_rounded[j] = round(fps[j], digits = 6)
    end
    fps_all = unique(sort(fps_rounded))

    # for certain v0 values and extra fp is detect. This is not correct
    if length(fps_all) > 3
        fps_all = fps_all[1:3]
    end

    # get left fixed point only
    fps_v = fps_all

    # return fps_v, fps_w, epsilon_crit, fps
    return fps_v

end
