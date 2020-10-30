# Solving ODEs using Forward Euler method

#Functions
w_inf1(v::Float64) = 10 ./ (1 + exp(-5 * (v - (1/3))))

# ODE's
dv(v::Float64,w::Float64,I_app::Float64) = (dt) .* (v -v^3/3 -w^2 +I_app)
dw(v::Float64,w::Float64,v0::Float64,w0::Float64,epsilon::Float64) = (dt) .* (epsilon*(w_inf1(v - v0) + (w0 - w)))


# Simulation using Euler
function FM_forward_euler_winf_new2(v0::Real, w0::Real, I_app::Real, epsilon::Real, v_init::Real, w_init::Real, t::Array{Float64,1})

    # Initial values
    v::Float64 = v_init
    vprev::Float64 = copy(v)
    w::Float64 = w_init

    # Pre-allocate
    vv = zeros(length(t))         # Vector in which V timecourse will be saved
    ww = zeros(length(t))         # Vector in which n timecourse will be saved

    # Euler method to find the next voltage value
    for i = 1:length(t)-1

        v += dv(vprev,w,I_app)              # Updates V
        w += dw(vprev,w,v0,w0,epsilon)      # Updates n using the current value of V, not the updated one

        vprev = copy(v)
        vv[i] = copy(v)             # V is saved in the vector VV
        ww[i] = copy(w)             # n is saved in the vector nn

    end

    return vv, ww
end


# Simulation with spike output
function FM_simulation_spkt_winf_new2(v0::Real, w0::Real, I_app::Real, epsilon::Real, v_init::Real, w_init::Real, t::Real)

    # initial values
    v::Float64=-1.          # Membrane potential
    vprev::Float64=-1.      # Previous value of membrane potential to update gating variables synchronously
    w::Float64=0.           # Potassium current activation variables

    spk_time = zeros(length(t))
    l=1

    # Update of variables at each time step
    for i = 1:length(t)

        v += dv(v,w,I_app)           # Updates V
        w += dw(vprev,w,v0,w0)      # Updates n using the current value of V, not the updated one

        #Get spikes
        if vprev < 0 && v >= 0
            spk_time[l] = t[i]
            l += 1
        end
        vprev = copy(V)
    end
    return spk_time
end
