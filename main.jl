include("circles_in_box_linear_spring_hysteresis.jl")
include("helpers.jl")

println(@__FILE__)
# domain description (box dimensions)
L = 1
H = 1
T = 0.15

# run settings
dt = 0.00025

# random initial conditions. if not given, random anyway and N = 25 by default
N = 25
initial_position, R = generate_circles(N, L, H)
initial_velocities = randn(Float64, (N, 2))
initial_accelerations = zeros(N, 2)
initial_angles = 2*pi*rand(Float64, N)
initial_omegas = 2*(2*pi*rand(Float64, N))
initial_angular_accelerations = zeros(N) 
initial_condition = (N, R, initial_position,
		     initial_velocities, initial_accelerations,
		     initial_angles, initial_omegas,
		     initial_angular_accelerations)

save_gif = false
verbose = false
# material is specified inside function
circles_in_box_linear_spring_hysteresis(L, H, T, dt, save_gif, verbose, initial_condition)
