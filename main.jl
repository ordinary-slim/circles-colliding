include("ball_in_box.jl")
println(@__FILE__)
# domain description
L = 1
H = 1

# run settings
T = 1
dt = 0.001

# single particle simulation

#initial_position = R .+ (1-2*R)*rand(1, 2)
#initial_velocity = rand(1, 2)
N = 150
save_gif = false
verbose = false
many_circles_in_box(N, L, H, T, dt, save_gif, verbose)
