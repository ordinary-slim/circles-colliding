using Plots
using LinearAlgebra
include("helpers.jl")
gr()
function circles_in_box_linear_spring_dashpot(material, L, H, T, dt, save_gif, verbose, initial_condition...)
	" linear spring with dashpot in normal direction"
	" linear spring in tangent direction"
	
	#unpack dictionary of material properties
	E = material["E"]# N/mm^2
	restitution_coeff = material["restitution_coeff"]
	mu_s = material["mu_s"]
	mu_roll = material["mu_roll"]
	tangential_stiff_ratio = material["tangential_stiff_ratio"]
	rho = material["rho"]

	epsilon =  1e-5# for comparison with small numbers

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
	angles = Array{Array{Float64, 1}, 1}(undef, N)
	angular_velocities = Array{Array{Float64, 1}, 1}(undef, N)
	angular_accelerations = Array{Array{Float64, 1}, 1}(undef, N)


	translational_kinetic_energy = sum(M.*sum(initial_velocities.^2, dims=2))/2
	rotational_kinetic_energy = sum(rot_intertia.*(initial_omegas.^2))/2
	println("Translational kinetic energy : ", 1000*translational_kinetic_energy, " mJ")
	println("Rotational kinetic energy : ", 1000*rotational_kinetic_energy, " mJ")
	println("Total kinetic energy of system : ", 1000*(translational_kinetic_energy+rotational_kinetic_energy), " mJ")
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
		angles[i] = zeros(number_iterations+1)
		angular_velocities[i] = zeros(number_iterations+1)
		angular_accelerations[i] = zeros(number_iterations+1)

		positions[i][1, :] = initial_positions[i, :]
		velocities[i][1, :] = initial_velocities[i, :]
		accelerations[i][1, :] = initial_accelerations[i, :]
		angles[i][1] = initial_angles[i]
		angular_velocities[i][1] = initial_omegas[i]
		angular_accelerations[i][1] = initial_angular_accelerations[i]
						       
		endpoints_x[2*(i+1)-1] = initial_positions[i, 1] - R[i]
		endpoints_x[2*(i+1)] = initial_positions[i, 1] + R[i]
		endpoints_y[2*(i+1)-1] = initial_positions[i, 2] - R[i]
		endpoints_y[2*(i+1)] = initial_positions[i, 2] + R[i]
	end


	# Visualization settings
	circle_resolution = 100
	theta = LinRange(0, 2*pi, circle_resolution)
	function ball_plot!(x, y, R, alpha, L, H)
		tmp = plot!(x.+ R*cos.(theta), y.+ R*sin.(theta), xlims = (0, L), ylims = (0, H), aspect_ratio =1, seriestype = [:shape], lw = 0, legend = false, framestyle = :none, grid = false, ticks = false, windowsize = (1200, 900), background_color_inside = :white, background_color_outside = :lightskyblue4, fillalpha = 0.6, linealpha = 0)
		return plot!([x+0.4*R*cos(alpha); x+0.9*R*cos(alpha)], [y+0.4*R*sin(alpha); y+0.9*R*sin(alpha)], width = 4, linealpha = 0.4)
	end
	my_plot = plot()
	for i=1:N
		if i==1
			my_plot = plot()
		end
		ball_plot!(positions[i][1, 1], positions[i][1, 2], R[i], angles[i][1], L, H)
	end
	display(my_plot)
	sleep(0.025)


	t = 0
	collision_prev_iter = Dict{Tuple{Int64, Int64}, Tuple{Array{Float64, 1}, Float64, Float64}}()
	#keys are currently colliding pairs, entries are position in the force-displacement plane (delta, F)
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
			x_1 = positions[i][iteration_number, :]
			x_2 = positions[j][iteration_number, :]
			v_1 = velocities[i][iteration_number, :]
			v_2 = velocities[j][iteration_number, :]
			omega_1 = angular_velocities[i][iteration_number]
			omega_2 = angular_velocities[j][iteration_number]
			overlap = R[i] + R[j] - norm(x_1 - x_2) #check if overlapped
			rel_vel = v_2 - v_1
			were_colliding = haskey(collision_prev_iter,  (i, j))
			just_occurred = false
			if (overlap>=0)
				normal = (x_2 - x_1)/norm(x_2 - x_1)#vector along centroid connecting line
				normal = reshape(normal, (2, ))# convert to column vector
				tangent = [-normal[2]; normal[1]]
				Q = [normal tangent]#from local to global

				if !were_colliding
					if verbose
						println("New collision between balls ", i, " and ", j, ".")
					end
					collision_prev_iter[(i, j)] = (normal, 0, 0)
					just_occurred = true
				end

				relative_pos = (x_2 - x_1)
				if just_occurred
					relative_pos_tm1 = relative_pos
					dTheta_1 = 0
					dTheta_2 = 0
				else
					relative_pos_tm1 = positions[j][iteration_number-1, :] - positions[i][iteration_number-1, :]
					dTheta_1 = (angles[i][iteration_number] - angles[i][iteration_number-1])
					dTheta_2 = (angles[j][iteration_number] - angles[j][iteration_number-1])
				end
				relative_displacement = relative_pos - relative_pos_tm1
				relative_displacement = Q'*relative_displacement#switch to local coordinates
				rel_vel = Q'*rel_vel
				overlap_dotted = -rel_vel[1]
				# right now  we assume that R-overlap/2>0
				distance2contactpoint_1 = (R[i]-overlap/2)
				distance2contactpoint_2 = (R[j]-overlap/2)
				relative_displacement[2] += -dTheta_2*distance2contactpoint_2 - dTheta_1*distance2contactpoint_1# add rolling to tangential displacement

				# deformation history
				normal_tm1, Fn_tm1, Ft_tm1 = collision_prev_iter[(i, j)]# t minus 1

				Fn, Ft = 0.0, 0.0
				# linear spring normal direction
			        k_nli, k_nlj = E*2*R[i], E*2*R[j]# taking as particle size the diameter
				k_nl = 1/(1/k_nli + 1/k_nlj)
				Fn += k_nl*overlap
				# dashpot normal direction
				effective_mass = 1/(1/M[i] + 1/M[j])
				c = sqrt(4*effective_mass*k_nl/(1 + (pi/log(restitution_coeff))^2))# EDEM version
				Fn += c*overlap_dotted

				# tangential forces
				# linear spring along tangential direction
				tangent_tm1 = [-normal_tm1[2]; normal_tm1[1]]
				k_t = tangential_stiff_ratio*k_nl
				pseudo_Ft = Ft_tm1*tangent_tm1 - k_t*relative_displacement[2]*tangent
				# Coulomb limit
				if (norm(pseudo_Ft)>mu_s*Fn)&&(norm(pseudo_Ft)>epsilon)
					pseudo_Ft = mu_s*Fn*pseudo_Ft/norm(pseudo_Ft)
				end
					
				Fn += dot(pseudo_Ft, normal)
				Ft += dot(pseudo_Ft, tangent)

				# torque from tangential forces
				Tao_1 = -distance2contactpoint_1*Ft
				Tao_2 = -distance2contactpoint_2*Ft
				# torque from rolling friction
				Tao_1 += -mu_roll*(R[i]-overlap/2)*Fn*sign(omega_1)
				Tao_2 += -mu_roll*(R[j]-overlap/2)*Fn*sign(omega_2)

				# update to current state
				collision_prev_iter[(i, j)] = (normal, Fn, Ft)

				a_1n = -Fn/M[i]
				a_2n = Fn/M[j]
				a_1t = -Ft/M[i]
				a_2t = Ft/M[j]
				angular_accelerations[i][iteration_number] += Tao_1/rot_intertia[i]
				angular_accelerations[j][iteration_number] += Tao_2/rot_intertia[j]
				accelerations[i][iteration_number, :] += a_1n*normal + a_1t*tangent
				accelerations[j][iteration_number, :] += a_2n*normal + a_2t*tangent
			 else
				if were_colliding
					if verbose
						println("Balls ", i, " and ", j, " are no longer in contact.")
					end
					delete!(collision_prev_iter, (i, j))
				end
			 end
		end
		####################################################################
		# update motion and bounding boxes
		t += dt
		#if verbose
			#println("Time :", t, " seconds.")
		#end
		for i=1:N
			#translational motion
			velocities[i][iteration_number+1, 1] = velocities[i][iteration_number, 1] + accelerations[i][iteration_number, 1]*dt
			velocities[i][iteration_number+1, 2] = velocities[i][iteration_number, 2] + accelerations[i][iteration_number, 2]*dt
			#positions[i][iteration_number+1, 1] = positions[i][iteration_number, 1] + velocities[i][iteration_number, 1]*dt
			#positions[i][iteration_number+1, 2] = positions[i][iteration_number, 2] + velocities[i][iteration_number, 2]*dt
			positions[i][iteration_number+1, 1] = positions[i][iteration_number, 1] + velocities[i][iteration_number, 1]*dt + accelerations[i][iteration_number, 1]*(dt^2)/2
			positions[i][iteration_number+1, 2] = positions[i][iteration_number, 2] + velocities[i][iteration_number, 2]*dt + accelerations[i][iteration_number, 2]*(dt^2)/2

			#rotational motion
			angular_velocities[i][iteration_number+1] = angular_velocities[i][iteration_number] + angular_accelerations[i][iteration_number]*dt
			angles[i][iteration_number+1] = angles[i][iteration_number] + angular_velocities[i][iteration_number]*dt + angular_accelerations[i][iteration_number]*(dt^2)/2
			if i==1
				my_plot = plot()
			end
			endpoints_x[2*(i+1)-1] = positions[i][iteration_number+1, 1] - R[i]
			endpoints_x[2*(i+1)] = positions[i][iteration_number+1, 1] + R[i]
			endpoints_y[2*(i+1)-1] = positions[i][iteration_number+1, 2] - R[i]
			endpoints_y[2*(i+1)] = positions[i][iteration_number+1, 2] + R[i]
			ball_plot!(positions[i][iteration_number+1, 1], positions[i][iteration_number+1, 2], R[i], angles[i][iteration_number + 1], L, H)
		end
		display(my_plot)
		sleep(0.025)
	end

	translational_kinetic_energy = 0
	rotational_kinetic_energy = 0
	for ith = 1:N
		translational_kinetic_energy += M[ith]*sum(velocities[ith][end, :].^2)/2
		rotational_kinetic_energy += rot_intertia[ith]*angular_velocities[ith][end]^2/2
	end
	
	println("Translational kinetic energy : ", 1000*translational_kinetic_energy, " mJ")
	println("Rotational kinetic energy : ", 1000*rotational_kinetic_energy, " mJ")
	println("Final kinetic energy of system : ", 1000*(translational_kinetic_energy + rotational_kinetic_energy), " mJ")

	if save_gif
		anim = @animate for iteration_number=1:number_iterations
			for i=1:N
				if i==1
					my_plot = plot()
				end
				ball_plot!(positions[i][iteration_number+1, 1], positions[i][iteration_number+1, 2], R[i], angles[i][iteration_number + 1], L, H)
			end
		end
		frames_per_second = 72
		name_o_file = string(N, "LinearSpringDashpot", frames_per_second, "fps.gif")
		gif(anim, name_o_file, fps = frames_per_second)
	end
end
