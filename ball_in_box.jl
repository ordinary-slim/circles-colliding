using Plots
using LinearAlgebra
<<<<<<< HEAD
include("helpers.jl")
gr()
function circles_in_box_linear_spring_hysteresis(N, L, H, T, dt, save_gif, verbose)
	" linear spring hysteresis contact model"
	E = 1000#N/mm^2
	restitution_coeff = 0.5
	mu_s = 0.4
	mu_roll = 0
	tangential_stiff_ratio = 1
	rho = 1

=======
gr()
function single_point_in_box(initial_position, initial_velocity, L, H, T, dt)
	number_iterations = Int(T/dt + 1)
	positions = Array{Float64, 2}(undef, (number_iterations, 2))
	velocities = Array{Float64, 2}(undef, (number_iterations, 2))

	positions[1, :] = initial_position
	velocities[1, :] = initial_velocity

	t = 0
	i = 0
	my_plot = scatter((positions[1, 1], positions[1, 2]), xlims = (0, L), ylims = (0, H))
	display(my_plot)
	sleep(0.025)
	while(t<=T)
		i+= 1
		x = positions[i, 1]
		y = positions[i, 2]
		v_x = velocities[i, 1]
		v_y = velocities[i, 2]

		## check for collision
		if (x<=0) | (x>= 1) # side wall collision
			v_x = -v_x
		end
		if (y<=0) | (y>= 1) # roof/ground collision
			v_y = -v_y
		end

		# update position
		x = x + v_x*dt
		y = y + v_y*dt
		if(t+dt>T)
			return
		else
			# update
			t += dt
			positions[i+1, 1] = x
			positions[i+1, 2] = y
			velocities[i+1, 1] = v_x
			velocities[i+1, 2] = v_y
			my_plot = scatter((positions[i, 1], positions[i, 2]), xlims = (0, L), ylims = (0, H))
			display(my_plot)
			sleep(0.025)
		end
	end
end
function many_points_in_box_no_interaction(N, L, H, T, dt)
	number_iterations = Int(T/dt + 1)
	positions = Array{Array{Float64, 2}, 1}(undef, N)
	velocities = Array{Array{Float64, 2}, 1}(undef, N)
	# initialize
	for i=1:N
		positions[i] = zeros(number_iterations, 2)
		velocities[i] = zeros(number_iterations, 2)
		positions[i][1, :] = rand(1, 2)
		velocities[i][1, :] = rand(1, 2)
	end

	t = 0
	my_plot = plot()
	plot!(xlims = (0, L), ylims = (0, H))
	for i=1:N
		my_plot = scatter!((positions[i][1, 1], positions[i][1, 2]), show=false, legend=false)
	end
	display(my_plot)
	sleep(0.025)
	my_plot = plot(xlims = (0, L), ylims = (0, H))
	j = 0
	while(t<=T)
		j+= 1
		for i=1:N
			x = positions[i][j, 1]
			y = positions[i][j, 2]
			v_x = velocities[i][j, 1]
			v_y = velocities[i][j, 2]

			## check for collision
			if (x<=0) | (x>= 1) # side wall collision
				v_x = -v_x
			end
			if (y<=0) | (y>= 1) # roof/ground collision
				v_y = -v_y
			end

			# update position
			x = x + v_x*dt
			y = y + v_y*dt
			if(t+dt>T)
				return
			else
				# update
				positions[i][j+1, 1] = x
				positions[i][j+1, 2] = y
				velocities[i][j+1, 1] = v_x
				velocities[i][j+1, 2] = v_y
				scatter!((positions[i][j, 1], positions[i][j, 2]), show=false, legend=false)
			end
		end
		t += dt
		display(my_plot)
		sleep(0.025)
		my_plot = plot(xlims = (0, L), ylims = (0, H))
	end
end
function single_circle_in_box(R, initial_position, initial_velocity, L, H, T, dt)
	" this time with gravity "
	# check if IC are out of bounds
	if (initial_position[1]<R || (L-initial_position[1])<R || initial_position[2]<R || (H-initial_position[2])<R)
		print("! The ball is already stuck in the wall!")
		return
	end
	number_iterations = Int(T/dt + 1)
	# Visualization settings
	circle_resolution = 201
	theta = LinRange(0, 2*pi, circle_resolution)
	circle_x, circle_y = R*cos.(theta), R*sin.(theta)
	function ball_plot(x, y)
		return plot(x, y, xlims = (0, L), ylims = (0, H), aspect_ratio =1, seriestype = [:shape], lw = 6, legend = false, framestyle = :box, grid = false, ticks = false, windowsize = (1200, 900), background_color_inside = :papayawhip, background_color_outside = :tomato4, fillcolor = :firebrick1, fillalpha = 0.9)
	end


	positions = Array{Float64, 2}(undef, (number_iterations, 2)) # center of circle
	velocities = Array{Float64, 2}(undef, (number_iterations, 2))

	positions[1, :] = initial_position
	velocities[1, :] = initial_velocity

	t = 0
	i = 0
	my_plot = ball_plot(positions[1, 1] .+ circle_x, positions[1, 2] .+ circle_y) 
	display(my_plot)
	sleep(0.025)
	while(t<=T)
		i+= 1
		x = positions[i, 1]
		y = positions[i, 2]
		v_x = velocities[i, 1]
		v_y = velocities[i, 2]
		a_x = 0
		a_y = -0.5

		# check for collision
		if (x<R || (L-x)<R) # side wall collision
			a_x = -2*v_x/dt # perfectly elastic shock, impulse calculation
		end
		if (y<R || (H-y)<R) # roof/ground collision
			a_y = -2*v_y/dt # perfectly elastic shock, impulse calculation
		end

		# update motion
		v_x = v_x + a_x*dt
		v_y = v_y + a_y*dt
		x = x + v_x*dt + a_x*(dt^2)/2
		y = y + v_y*dt + a_y*(dt^2)/2
		if(t+dt>T)
			return
		else
			# update
			t += dt
			positions[i+1, 1] = x
			positions[i+1, 2] = y
			velocities[i+1, 1] = v_x
			velocities[i+1, 2] = v_y
			my_plot = ball_plot(positions[i, 1] .+ circle_x, positions[i, 2] .+ circle_y) 
			display(my_plot)
			sleep(0.025)
		end
	end
end
function generate_circles(N, L, H, rho)
	" generate non-overlapping circles in box of width L and height H "
	" boxes are of homogeneneous material of density rho "
	# upper and lower bounds for the radius
	upper_bound = 1
	lower_bound = sqrt(L*H/(2*N))/4
	radii = lower_bound .+ (upper_bound-lower_bound)*rand(Float64, N)
	centers = rand(Float64, (N, 2))
	centers[:, 1] *= L
	centers[:, 2] *= H
	masses = zeros(N)
	for i=1:N
		if i==1 # only do wall check
			while ((centers[i, 1] <= radii[i])||((L -centers[i, 1]) <= radii[i])||(centers[i, 2] <= radii[i])||((H -centers[i, 2]) <= radii[i]))
				radii[i] = lower_bound + upper_bound*rand()
				centers[i, :] = rand(Float64, (1, 2))
				centers[i, 1] *= L
				centers[i, 2] *= H
			end
		else
			ball_approved = false
			while ~ball_approved
				for j=1:(i-1)
					if (norm(centers[j, :] - centers[i, :]) <= (radii[j] + radii[i])||((centers[i, 1] <= radii[i])||((L -centers[i, 1]) <= radii[i])||(centers[i, 2] <= radii[i])||((H -centers[i, 2]) <= radii[i]))) #test failed
						radii[i] = lower_bound + upper_bound*rand()
						centers[i, :] = rand(Float64, (1, 2))
						centers[i, 1] *= L
						centers[i, 2] *= H
						break
					elseif j==(i-1) # last test validated
						ball_approved = true
					end
				end
			end
		end
		masses[i] = rho*pi*radii[i]^2
	end
	return centers, radii, masses
end

function overlapping_circles_check(x_1, x_2, r_1, r_2)
	" circles "
	return ((norm(x_1 - x_2) <= r_1 + r_2)) #check if overlapped
end

function getting_closer_check(x_1, x_2, v_1, v_2)
	" general "
	return (dot(v_2 - v_1, x_2 - x_1) < 0 ) #check if getting closer
end

function sweep_n_prune(endpoints, indices)
	"
	endpoints : 1D array of Floats that contains the endpoints of the objects' bounding boxes in unsorted order
	indices : 1D array of Ints. indices[i] = index of object who endpoints[i] belongs to
	sweep_n_prune is a BPC algorithm that detects overlaps between the bounding boxes
	left and right wall correspond to indices 0 and N+1, where N is the number of objects
	"
	M = Int64(length(indices)/2) # M = N+2
	# first sweep
	p = sortperm(endpoints, alg=InsertionSort)
	# second sweep
	PCS = Set{Tuple{Int64, Int64}}()
	AS = Set{Int64}()
	LWCS = Set{Int64}()
	RWCS = Set{Int64}()
	is_active = zeros(Int8, M)
	for object in indices[p]
		if is_active[object+1] == 0
			is_active[object+1] = 1
			# for all already active items, possible collision
			for active_object in AS
				i, j = sort([object, active_object])
				if i == 0
					push!(LWCS, j)
				elseif j == M-1
					push!(RWCS, i)
				else
					push!(PCS, (i, j))
				end
			end
			push!(AS, object)
		else
			is_active[object+1] = 0
			setdiff!(AS, object)
		end
	end
	return p, LWCS, RWCS, PCS
end

function many_circles_in_box(N, L, H, T, dt, save_gif, verbose)
	" first implementation with check of all balls for collisions "
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
	number_iterations = Int64(div(T, dt, RoundUp))
	positions = Array{Array{Float64, 2}, 1}(undef, N)
	velocities = Array{Array{Float64, 2}, 1}(undef, N)
	accelerations = Array{Array{Float64, 2}, 1}(undef, N)
<<<<<<< HEAD
	angles = Array{Array{Float64, 1}, 1}(undef, N)
	angular_velocities = Array{Array{Float64, 1}, 1}(undef, N)
	angular_accelerations = Array{Array{Float64, 1}, 1}(undef, N)

	# initialize
	initial_position, R, M = generate_circles(N, L, H, rho)
	rot_intertia = M.*(R.^2)/2# moments of intertia
	initial_velocity = randn(Float64, (N, 2))
	initial_angles = 2*pi*rand(Float64, N)
	initial_angular_velocities = 2*pi*rand(Float64, N)

	# shearing test
	#N = 2
	#initial_position = [0.4 0.35; 0.6 0.65]
	#R = [0.15 0.15]
	#M = pi*R.^2*rho
	#rot_intertia = M.*(R.^2)/2# moments of intertia
	#initial_velocity = [0 1.0; 0 -1.0]
	#initial_angles = [0.0 0.0]
	#initial_angular_velocities = [0.0 0.0]

	kinetic_energy = sum(M.*sum(initial_velocity.^2, dims=2))/2
	kinetic_energy += sum(rot_intertia.*(initial_angular_velocities.^2))/2
=======
	# initialize
	rho = 1
	initial_position, R, M = generate_circles(N, L, H, rho)
	initial_velocity = randn(Float64, (N, 2))
	kinetic_energy = sum(M.*(sum(initial_velocity.^2, dims=2).^0.5))
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
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
<<<<<<< HEAD
		angles[i] = zeros(number_iterations+1)
		angular_velocities[i] = zeros(number_iterations+1)
		angular_accelerations[i] = zeros(number_iterations+1)

		positions[i][1, :] = initial_position[i, :]
		velocities[i][1, :] = initial_velocity[i, :]
		angles[i][1] = initial_angles[i]
		angular_velocities[i][1] = initial_angular_velocities[i]
						       
=======
		positions[i][1, :] = initial_position[i, :]
		velocities[i][1, :] = initial_velocity[i, :]
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
		endpoints_x[2*(i+1)-1] = initial_position[i, 1] - R[i]
		endpoints_x[2*(i+1)] = initial_position[i, 1] + R[i]
		endpoints_y[2*(i+1)-1] = initial_position[i, 2] - R[i]
		endpoints_y[2*(i+1)] = initial_position[i, 2] + R[i]
	end


	# Visualization settings
	circle_resolution = 50
	theta = LinRange(0, 2*pi, circle_resolution)
<<<<<<< HEAD
	function ball_plot!(x, y, R, alpha, L, H)
		plot!(x.+ R*cos.(theta), y.+ R*sin.(theta), xlims = (0, L), ylims = (0, H), aspect_ratio =1, seriestype = [:shape], lw = 0, legend = false, framestyle = :none, grid = false, ticks = false, windowsize = (1200, 900), background_color_inside = :white, background_color_outside = :lightskyblue4, fillalpha = 0.6, linealpha = 0)
		return plot!([x; x+R*cos(alpha)], [y; y+R*sin(alpha)])
=======
	function ball_plot!(x, y, L, H)
		return plot!(x, y, xlims = (0, L), ylims = (0, H), aspect_ratio =1, seriestype = [:shape], lw = 0, legend = false, framestyle = :none, grid = false, ticks = false, windowsize = (1200, 900), background_color_inside = :white, background_color_outside = :lightskyblue4, fillalpha = 0.6, linealpha = 0)
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
	end
	my_plot = plot()
	for i=1:N
		if i==1
			my_plot = plot()
		end
<<<<<<< HEAD
		ball_plot!(positions[i][1, 1], positions[i][1, 2], R[i], angles[i][1], L, H)
=======
		ball_plot!(positions[i][1, 1] .+ R[i]*cos.(theta), positions[i][1, 2] .+ R[i]*sin.(theta), L, H)
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
	end
	display(my_plot)
	sleep(0.025)


	t = 0
<<<<<<< HEAD
	collision_prev_iter = Dict{Tuple{Int64, Int64}, Tuple{Array{Float64, 1}, Float64, Float64}}()
	#keys are currently colliding pairs, entries are position in the force-displacement plane (delta, F)
=======
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
	for iteration_number=1:number_iterations
		if(t+dt>T)
			dt = T - t
		end
		# calculate forces then update motion

		# check for collisions
		# collision relation is "anti-reflexive" and symmetric
<<<<<<< HEAD
		relative_p_x, LWCS, RWCS, PCS_x = sweep_n_prune(endpoints_x[p_x], indices_x[p_x])
		relative_p_y, BWCS, TWCS, PCS_y = sweep_n_prune(endpoints_y[p_y], indices_y[p_y])
		p_x = p_x[relative_p_x]
		p_y = p_y[relative_p_y]
		PCS = intersect(PCS_x, PCS_y)
=======
		println("Sweep n prune in x-axis:")
		@time relative_p_x, LWCS, RWCS, PCS_x = sweep_n_prune(endpoints_x[p_x], indices_x[p_x])
		println("Sweep n prune in y-axis:")
		@time relative_p_y, BWCS, TWCS, PCS_y = sweep_n_prune(endpoints_y[p_y], indices_y[p_y])
		p_x = p_x[relative_p_x]
		p_y = p_y[relative_p_y]
		println("Intersection:")
		@time PCS = intersect(PCS_x, PCS_y)
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
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
<<<<<<< HEAD
			x_1 = positions[i][iteration_number, :]
			x_2 = positions[j][iteration_number, :]
			v_1 = velocities[i][iteration_number, :]
			v_2 = velocities[j][iteration_number, :]
			omega_1 = angular_velocities[i][iteration_number]
			omega_2 = angular_velocities[j][iteration_number]
			overlap = R[i] + R[j] - norm(x_1 - x_2) #check if overlapped
			rel_vel = v_2 - v_1
			were_colliding = haskey(collision_prev_iter,  (i, j))
			if (overlap>=0)
				if !were_colliding
					if verbose
						println("New collision between balls ", i, " and ", j, ".")
					end
					collision_prev_iter[(i, j)] = ([0.0; 0.0], 0, 0)
				end
				normal = (x_2 - x_1)/norm(x_2 - x_1)#vector along centroid connecting line
				normal = reshape(normal, (2, ))# convert to column vector
				tangent = [-normal[2]; normal[1]]
				Q = [normal tangent]#from local to global

				relative_pos_tm1 = positions[j][iteration_number-1, :] - positions[i][iteration_number-1, :]
				relative_pos = (x_2 - x_1)
				relative_displacement = relative_pos - relative_pos_tm1
				relative_displacement = Q'*relative_displacement#switch to local coordinates

				# normal forces
				# Linear spring with hysteresis along normal direction
				normal_tm1, Fn_tm1, Ft_tm1 = collision_prev_iter[(i, j)]# t minus 1

			        k_nli, k_nlj = E*R[i], E*R[j]
				k_nl = 1/(1/k_nli + 1/k_nlj)
				k_nul = k_nl/restitution_coeff^2

				loading = (dot(rel_vel, normal)<=0)
				if loading
					Fn = min(k_nl*overlap, Fn_tm1 - relative_displacement[1]*k_nul)
				else
					Fn = max(Fn_tm1 - relative_displacement[1]*k_nul, 0)
				end

				# tangential forces
				# linear spring along tangential direction
				tangent_tm1 = [-normal_tm1[2]; normal_tm1[1]]
				k_t = tangential_stiff_ratio*k_nl
				pseudo_Ft = Ft_tm1*tangent_tm1 - k_t*relative_displacement[2]*tangent
				pseudo_Ft = max(norm(pseudo_Ft), mu_s*Fn)*pseudo_Ft/norm(pseudo_Ft)
				Fn += dot(pseudo_Ft, normal)
				Ft = dot(pseudo_Ft, tangent)

				# torque from tangential forces
				Tao_1 = -(R[i]-overlap/2)/2*Ft
				Tao_2 = -(R[j]-overlap/2)/2*Ft
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
=======
			overlapped = overlapping_circles_check(positions[i][iteration_number, :], positions[j][iteration_number, :], R[i], R[j])
			getting_closer = getting_closer_check(positions[i][iteration_number, :], positions[j][iteration_number, :], velocities[i][iteration_number, :], velocities[j][iteration_number, :])
			if overlapped&&getting_closer
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
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
		end
		####################################################################
		# update motion and bounding boxes
		t += dt
<<<<<<< HEAD
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
=======
		println("Time :", t, " seconds.")
		for i=1:N
			velocities[i][iteration_number+1, 1] = velocities[i][iteration_number, 1] + accelerations[i][iteration_number, 1]*dt
			velocities[i][iteration_number+1, 2] = velocities[i][iteration_number, 2] + accelerations[i][iteration_number, 2]*dt
			positions[i][iteration_number+1, 1] = positions[i][iteration_number, 1] + velocities[i][iteration_number, 1]*dt + accelerations[i][iteration_number, 1]*(dt^2)/2
			positions[i][iteration_number+1, 2] = positions[i][iteration_number, 2] + velocities[i][iteration_number, 2]*dt + accelerations[i][iteration_number, 2]*(dt^2)/2
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
			if i==1
				my_plot = plot()
			end
			endpoints_x[2*(i+1)-1] = positions[i][iteration_number+1, 1] - R[i]
			endpoints_x[2*(i+1)] = positions[i][iteration_number+1, 1] + R[i]
			endpoints_y[2*(i+1)-1] = positions[i][iteration_number+1, 2] - R[i]
			endpoints_y[2*(i+1)] = positions[i][iteration_number+1, 2] + R[i]
<<<<<<< HEAD
			ball_plot!(positions[i][iteration_number+1, 1], positions[i][iteration_number+1, 2], R[i], angles[i][iteration_number + 1], L, H)
=======
			ball_plot!(positions[i][iteration_number+1, 1] .+ R[i]*cos.(theta), positions[i][iteration_number+1, 2] .+ R[i]*sin.(theta), L, H)
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
		end
		display(my_plot)
		sleep(0.025)
	end

<<<<<<< HEAD
	kinetic_energy = 0
	for ith = 1:N
		kinetic_energy += M[ith]*sum(velocities[ith][end, :].^2)/2
		kinetic_energy += rot_intertia[ith]*angular_velocities[ith][end]^2/2
	end
	println("Final kinetic energy of system : ", 1000*kinetic_energy, " mJ")

=======
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
	if save_gif
		anim = @animate for iteration_number=1:number_iterations
			for i=1:N
				if i==1
					my_plot = plot()
				end
<<<<<<< HEAD
				ball_plot!(positions[i][iteration_number+1, 1], positions[i][iteration_number+1, 2], R[i], angles[i][iteration_number + 1], L, H)
			end
		end
		gif(anim, "anim_fps15.gif", fps = 24)
=======
				ball_plot!(positions[i][iteration_number, 1] .+ R[i]*cos.(theta), positions[i][iteration_number, 2] .+ R[i]*sin.(theta), L, H)
			end
		end
		gif(anim, "anim_fps15.gif", fps = 15)
>>>>>>> cac67b37dc08ec3dbf9110fac25ead1f46b6eb48
	end
end
