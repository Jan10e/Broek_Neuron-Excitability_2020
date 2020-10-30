# Packages
using LinearAlgebra

# Functions
include("FM_forward_euler_step.jl")
include("FM_fixed_points_multi.jl")

# Computes spiking frequency
function frequency(y::Array{Float64,1})
    spk::Int64 = 0
    N0::Int64 = 0
    N1::Int64 = 0
    freq::Float64 = 0
    ISI::Float64 = 0

    for i = 1:length(y)-1
        if spk == 0 && y[i] < 0 && y[i+1] >= 0
            spk += 1
            N0 = copy(i)
        elseif spk == 1 && y[i] < 0 && y[i+1] >= 0
            N1 = copy(i)
            break
        end
    end
    if N1 > 0
        ISI = (N1 - N0) * dt
        freq = 1000 / ISI
    end
    return freq
end

# Check stability right fixed point
function check_rfp(v0::Real, w0::Real, I_app::Real, epsilon::Real)
    winf1(v::Real) = 10 / (1 + exp(-5 * (v - (1 / 3))))
    winf3(v::Float64, v0::Float64, w0::Float64) =
        10 / (1 + exp(-5 * (v - (1 / 3) - v0))) + w0
    steady_state(x) = (x[1] - (1 / 3) * x[1]^3 - (winf1(x[1] - v0) + w0)^2 + I_app)
    jacobian(fps_v::Real, fps_w::Real) = [
        [1 - fps_v^2, -2 * fps_w]
        [
            (epsilon * 50 * exp(5 * v0 - 5 * fps_v + 5 / 3)) /
            (exp(5 * v0 - 5 * fps_v + 5 / 3) + 1)^2,
            -epsilon,
        ]
    ]

    stable = nothing

    fps_v = FM_fixed_points_multi(steady_state)
    fp_v = fps_v[end] #get the rfp
    fp_w = winf3(fp_v, v0, w0)

    lambda = Array{Real,1}(undef, 2)
    jac = jacobian(fp_v, fp_w)
    jac_matrix = [jac[1] jac[2]; jac[3] jac[4]]
    lambda = eigvals(jac_matrix)

    if real(lambda[1]) < 0 && real(lambda[2]) < 0
        stable = 1
    else
        stable = 0
    end
    return stable
end


#Computes the fIcurve
function ficurve_step_down(
    v0::Real,
    w0::Real,
    epsilon::Real,
    I_app::Real,
    I_step_range::Vector{Float64},
)
    f = zeros(length(I_step_range))
    for i = 1:length(I_step_range)
        println(I_app - I_step_range[i])
        stable_rfp = check_rfp(v0, w0, -I_step_range[i], epsilon)
        if stable_rfp == 0
            (v, w) = FM_forward_euler_step(
                v0,
                w0,
                epsilon,
                I_app,
                T_step_start,
                T_step_stop,
                -I_step_range[i],
            )
            global v = v
            f[i] = frequency(v[Ttrans:end])
        elseif stable_rfp == 1
            f[i] = NaN
        end
    end
    return f, v
end
