include("ball_in_box.jl")
# domain description
L = 1
H = 1

# run settings
T = 2.5
dt = 0.005

# single particle simulation

#initial_position = R .+ (1-2*R)*rand(1, 2)
#initial_velocity = rand(1, 2)
N = 50
save_gif = false
verbose = false
many_circles_in_box(N, L, H, T, dt, save_gif, verbose)
