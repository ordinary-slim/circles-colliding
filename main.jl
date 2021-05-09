include("helpers.jl")
include("circles_in_box_perfectly_elastic.jl")
include("circles_in_box_linear_spring_hysteresis.jl")
include("circles_in_box_linear_spring_dashpot.jl")
include("circles_in_box_hertzian_spring_dashpot.jl")

println(@__FILE__)
# domain description (box dimensions)
L = 1
H = 1
T = 0.2

# material properties
material = Dict{String, Float64}()# note that without the "()", material is the constructor
material["rho"] = 1
material["E"] = 10000# N/mm^2
material["mu_s"] = 0.5
material["mu_roll"] = 0
material["tangential_stiff_ratio"] = 1
material["restitution_coeff"] = 0.8
material["poisson"] = 0.25

# run settings
dt = 0.0001

# random initial conditions. if not given, random anyway and N = 25 by default
N = 30
initial_position, R = generate_circles(N, L, H)
initial_velocities = 4*randn(Float64, (N, 2))
initial_accelerations = zeros(N, 2)
initial_angles = 2*pi*rand(Float64, N)
initial_omegas = 24*(2*pi*rand(Float64, N))
initial_angular_accelerations = zeros(N) 
initial_condition = (N, R, initial_position,
		     initial_velocities, initial_accelerations,
		     initial_angles, initial_omegas,
		     initial_angular_accelerations)

save_gif = true
verbose = false

#circles_in_box_perfectly_elastic(material, L, H, T, dt, save_gif, verbose, initial_condition)
#circles_in_box_linear_spring_hysteresis(material, L, H, T, dt, save_gif, verbose, initial_condition)
circles_in_box_linear_spring_dashpot(material, L, H, T, dt, save_gif, verbose, initial_condition)
#circles_in_box_hertzian_spring_dashpot(material, L, H, T, dt, save_gif, verbose, initial_condition)
