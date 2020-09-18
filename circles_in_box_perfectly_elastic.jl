using LinearAlgebra
using Plots
gr()
include("helpers.jl")
function circles_in_box_perfectly_elastic(material, L, H, T, dt, save_gif, verbose, initial_condition...)
	rho = material["rho"]
	if isempty(initial_condition)
		# initialize randomly
		N = 25
		initial_positions, R = generate_circles(N, L, H)
		M = pi*R.^2*rho
		rot_intertia = M.*(R.^2)/2# moments of intertia
		initial_velocities = randn(Float64, (N, 2))
		initial_accelerations = zeros(N, 2)
		initial_angles = 2*pi*rand(Float64, N)
		initial_omegas = 2*(2*pi*rand(Float64, N))
		initial_angular_accelerations = zeros(N) 

		# shearing test
		#N = 2
		#initial_positions = [0.4 0.35; 0.6 0.65]
		#R = [0.15 0.15]
		#M = pi*R.^2*rho
		#rot_intertia = M.*(R.^2)/2# moments of intertia
		#initial_velocities = [0 1.0; 0 -1.0]
		#initial_accelerations = zeros(N, 2)
		#initial_angles = [0.0 0.0]
		#initial_omegas = [-5.0 -5.0]
		#initial_angular_accelerations = zeros(N) 

		# rolling test
		#N = 2
		#initial_positions = [0.4 0.5; 0.6 0.5]
		#R = [0.15 0.15]
		#M = pi*R.^2*rho
		#rot_intertia = M.*(R.^2)/2# moments of intertia
		#initial_velocities = [0.0 0.0; 0.0 0.0]
		#initial_accelerations = zeros(N, 2)
		#initial_angles = [0.0 pi]
		#initial_omegas = [5.0 0.0]
		#initial_angular_accelerations = zeros(N) 
	else
		initial_condition = initial_condition[1]
		N = initial_condition[1]
		R = initial_condition[2]
		initial_positions, initial_velocities, initial_accelerations = initial_condition[3:5]
		initial_angles, initial_omegas, initial_angular_accelerations = initial_condition[6:8]

		M = pi*R.^2*rho
		rot_intertia = M.*(R.^2)/2# moments of intertia
	end

	number_iterations = Int64(div(T, dt, RoundUp))
	positions = Array{Array{Float64, 2}, 1}(undef, N)
	velocities = Array{Array{Float64, 2}, 1}(undef, N)
	accelerations = Array{Array{Float64, 2}, 1}(undef, N)

	kinetic_energy = sum(M.*sum(initial_velocities.^2, dims=2))/2
	println("Total kinetic energy of system : ", 1000*kinetic_energy, " mJ")
	# endpoints and indices of bounding boxes
	# indices 0 and N+1 correspond to left and right walls
	# the endpoints and indices of object k rest at the indices 2*(k+1)-1 and 2*(k+1)
	indices_x = repeat(0:N+1, inner=2)
	indices_y = repeat(0:N+1, inner=2)
	p_x = 1:2*(N+2)
	p_y = 1:2*(N+2)
	endpoints_x = zeros(2*N+4) # +4 for the two endpoints of ea wall
	endpoints_y = zeros(2*N+4)
	# wall endpoints
	endpoints_x[1:2] = [-Inf, 0]
	endpoints_x[2*N+3:2*N+4] = [L, +Inf]
	endpoints_y[1:2] = [-Inf, 0]
	endpoints_y[2*N+3:2*N+4] = [H, +Inf]
	# initialize motion and bounding boxes
	for i=1:N
		positions[i] = zeros(number_iterations+1, 2)
		velocities[i] = zeros(number_iterations+1, 2)
		accelerations[i] = zeros(number_iterations+1, 2)
		positions[i][1, :] = initial_positions[i, :]
		velocities[i][1, :] = initial_velocities[i, :]
		accelerations[i][1, :] = initial_accelerations[i, :]
		endpoints_x[2*(i+1)-1] = initial_positions[i, 1] - R[i]
		endpoints_x[2*(i+1)] = initial_positions[i, 1] + R[i]
		endpoints_y[2*(i+1)-1] = initial_positions[i, 2] - R[i]
		endpoints_y[2*(i+1)] = initial_positions[i, 2] + R[i]
	end


	# Visualization settings
	circle_resolution = 50
	theta = LinRange(0, 2*pi, circle_resolution)
	function ball_plot!(x, y, L, H)
		return plot!(x, y, xlims = (0, L), ylims = (0, H), aspect_ratio =1, seriestype = [:shape], lw = 0, legend = false, framestyle = :none, grid = false, ticks = false, windowsize = (1200, 900), background_color_inside = :white, background_color_outside = :lightskyblue4, fillalpha = 0.6, linealpha = 0)
	end
	my_plot = plot()
	for i=1:N
		if i==1
			my_plot = plot()
		end
		ball_plot!(positions[i][1, 1] .+ R[i]*cos.(theta), positions[i][1, 2] .+ R[i]*sin.(theta), L, H)
	end
	display(my_plot)
	sleep(0.025)


	t = 0
	for iteration_number=1:number_iterations
		if(t+dt>T)
			dt = T - t
		end
		# calculate forces then update motion

		# check for collisions
		# collision relation is "anti-reflexive" and symmetric
		relative_p_x, LWCS, RWCS, PCS_x = sweep_n_prune(endpoints_x[p_x], indices_x[p_x])
		relative_p_y, BWCS, TWCS, PCS_y = sweep_n_prune(endpoints_y[p_y], indices_y[p_y])
		p_x = p_x[relative_p_x]
		p_y = p_y[relative_p_y]
		PCS = intersect(PCS_x, PCS_y)
		for i in LWCS #left wall collision set
			if (velocities[i][iteration_number, 1] <0)
				if verbose
					print("Ball ", i, " collided with the left wall!\n")
				end
				accelerations[i][iteration_number, 1] += -2*velocities[i][iteration_number, 1]/dt # perfectly elastic shock, impulse calculation
			end
		end
		for i in RWCS #right wall collision set
			if (velocities[i][iteration_number, 1] >0)
				if verbose
					print("Ball ", i, " collided with the right wall!\n")
				end
				accelerations[i][iteration_number, 1] += -2*velocities[i][iteration_number, 1]/dt
			end
		end
		for i in BWCS #bottom wall collision set
			if (velocities[i][iteration_number, 2] <0)
				if verbose
					print("Ball ", i, " collided with the bottom wall!\n")
				end
				accelerations[i][iteration_number, 2] += -2*velocities[i][iteration_number, 2]/dt
			end
		end
		for i in TWCS #right wall collision set
			if (velocities[i][iteration_number, 2] >0)
				if verbose
					print("Ball ", i, " collided with the top wall!\n")
				end
				accelerations[i][iteration_number, 2] += -2*velocities[i][iteration_number, 2]/dt
			end
		end
		for pair in PCS
			i = pair[1]
			j = pair[2]
			overlap = R[i] + R[j] - norm(positions[i][iteration_number, :] - positions[j][iteration_number, :]) #check if overlapped
			getting_closer = getting_closer_check(positions[i][iteration_number, :], positions[j][iteration_number, :], velocities[i][iteration_number, :], velocities[j][iteration_number, :])
			if (overlap>=0)&&getting_closer
				if verbose
					print("Ball ", i, " and ball ", j, " just collided !\n")
				end
				# Hypothesis : forces only along normal to contact plane
				x_1 = positions[i][iteration_number, :]
				x_2 = positions[j][iteration_number, :]
				n = (x_2 - x_1)/norm(x_2 - x_1)
				v_1n = dot(velocities[i][iteration_number, :], n)
				v_2n = dot(velocities[j][iteration_number, :], n)
				a_1n = 2*M[j]/((M[i] + M[j])*dt)*(v_2n - v_1n)
				a_2n = 2*M[i]/((M[i] + M[j])*dt)*(v_1n - v_2n)
				accelerations[i][iteration_number, :] += a_1n*n
				accelerations[j][iteration_number, :] += a_2n*n
			end
		end
		####################################################################
		# update motion and bounding boxes
		t += dt
		if verbose
			println("Time :", t, " seconds.")
		end
		for i=1:N
			velocities[i][iteration_number+1, 1] = velocities[i][iteration_number, 1] + accelerations[i][iteration_number, 1]*dt
			velocities[i][iteration_number+1, 2] = velocities[i][iteration_number, 2] + accelerations[i][iteration_number, 2]*dt
			positions[i][iteration_number+1, 1] = positions[i][iteration_number, 1] + velocities[i][iteration_number, 1]*dt
			positions[i][iteration_number+1, 2] = positions[i][iteration_number, 2] + velocities[i][iteration_number, 2]*dt
			#positions[i][iteration_number+1, 1] = positions[i][iteration_number, 1] + velocities[i][iteration_number, 1]*dt + accelerations[i][iteration_number, 1]*(dt^2)/2
			#positions[i][iteration_number+1, 2] = positions[i][iteration_number, 2] + velocities[i][iteration_number, 2]*dt + accelerations[i][iteration_number, 2]*(dt^2)/2
			if i==1
				my_plot = plot()
			end
			endpoints_x[2*(i+1)-1] = positions[i][iteration_number+1, 1] - R[i]
			endpoints_x[2*(i+1)] = positions[i][iteration_number+1, 1] + R[i]
			endpoints_y[2*(i+1)-1] = positions[i][iteration_number+1, 2] - R[i]
			endpoints_y[2*(i+1)] = positions[i][iteration_number+1, 2] + R[i]
			ball_plot!(positions[i][iteration_number+1, 1] .+ R[i]*cos.(theta), positions[i][iteration_number+1, 2] .+ R[i]*sin.(theta), L, H)
		end
		display(my_plot)
		sleep(0.025)
	end

	kinetic_energy = 0
	for ith = 1:N
		kinetic_energy += M[ith]*sum(velocities[ith][end, :].^2)/2
	end
	println("Final kinetic energy of system : ", 1000*kinetic_energy, " mJ")

	if save_gif
		anim = @animate for iteration_number=1:number_iterations
			for i=1:N
				if i==1
					my_plot = plot()
				end
				ball_plot!(positions[i][iteration_number, 1] .+ R[i]*cos.(theta), positions[i][iteration_number, 2] .+ R[i]*sin.(theta), L, H)
			end
		end
		frames_per_second = 60
		name_o_file = string(N, "PerfectlyElastic", frames_per_second, "fps.gif")
		gif(anim, name_o_file, fps = 24)
	end
end
