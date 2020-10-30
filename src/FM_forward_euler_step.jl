#Functions
w_inf1(v::Float64) = 10 ./ (1 + exp(-5 * (v - (1/3))))

# ODE's
dv(v::Float64,w::Float64,I_app::Float64,I_app_step::Float64) = (dt) .* (v -v^3/3 -w^2 + I_app + I_app_step)
dw(v::Float64,w::Float64,v0::Float64,w0::Float64,epsilon::Float64) = (dt) .* (epsilon*(w_inf1(v - v0) + (w0 - w)))


# Simulation using Euler
function FM_forward_euler_step(v0::Real, w0::Real, epsilon::Real, I_app::Real, T_step_start::Float64, T_step_stop::Float64, I_step::Float64)

    # Initial values
    v::Float64 = v_init
    vprev::Float64 = copy(v)
    w::Float64 = w_init

    # step start & stop
    Tstart::Int64 = convert(Int64,T_step_start/dt)
    Tstop::Int64 = convert(Int64,T_step_stop/dt)

    # Pre-allocate
    vv = zeros(length(t))         # Vector in which V timecourse will be saved
    ww = zeros(length(t))         # Vector in which n timecourse will be saved

    # Euler method with step method
    for i = 1:Tdt
        if i >= Tstart && i <= Tstop
            I_app_step = I_step
        else
            I_app_step = 0.
        end

        v += dv(vprev, w, I_app, I_app_step)              # Updates V
        w += dw(vprev, w, v0, w0, epsilon)      # Updates n using the current value of V, not the updated one

        vprev = copy(v)
        vv[i] = copy(v)             # V is saved in the vector VV
        ww[i] = copy(w)             # n is saved in the vector nn

    end

    return vv, ww
end
