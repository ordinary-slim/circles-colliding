using Plots
using LinearAlgebra
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
function generate_circles(N, L, H)
	" generate non-overlapping circles in box of width L and height H "
	#upper_bound = sqrt(L*H/N)/4
	upper_bound = 1
	lower_bound = sqrt(L*H/(2*N))/4
	radii = lower_bound .+ upper_bound*rand(Float64, N)
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
		if i==1
			masses[i] = 1
		else
			masses[i] = radii[i]^2/radii[1]^2
		end
	end
	return centers, radii, masses
end

function overlapping_check(x_1, x_2, r_1, r_2)
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
	"
	N = Int64(length(indices)/2)
	# first sweep
	p = sortperm(endpoints, alg=InsertionSort)
	# second sweep
	PCS = Set{Tuple{Int64, Int64}}()
	AS = Set{Int64}()
	is_active = zeros(Int8, N)
	for object in indices[p]
		if is_active[object] == 0
			is_active[object] = 1
			# for all already active items, possible collision
			for active_object in AS
				if active_object<object
					push!(PCS, (active_object, object))
				else
					push!(PCS, (object, active_object))
				end
			end
			push!(AS, object)
		else
			is_active[object] = 0
			setdiff!(AS, object)
		end
	end
	return p, PCS
end

function many_circles_in_box(N, L, H, T, dt, save_gif, verbose)
	" first implementation with check of all balls for collisions "
	number_iterations = ceil(Int64, T/dt + 1)
	positions = Array{Array{Float64, 2}, 1}(undef, N)
	velocities = Array{Array{Float64, 2}, 1}(undef, N)
	accelerations = Array{Array{Float64, 2}, 1}(undef, N)
	# initialize
	@time initial_position, R, M = generate_circles(N, L, H)
	initial_velocity = randn(Float64, (N, 2))
	# endpoints and indices of bounding boxes
	indices_x = repeat(1:N, inner=2)
	indices_y = repeat(1:N, inner=2)
	p_x = 1:2*N
	p_y = 1:2*N
	endpoints_x = zeros(2*N)
	endpoints_y = zeros(2*N)
	for i=1:N
		positions[i] = zeros(number_iterations, 2)
		velocities[i] = zeros(number_iterations, 2)
		accelerations[i] = zeros(number_iterations, 2)
		positions[i][1, :] = initial_position[i, :]
		velocities[i][1, :] = initial_velocity[i, :]
		endpoints_x[2*i-1] = initial_position[i, 1] - R[i]
		endpoints_x[2*i] = initial_position[i, 1] + R[i]
		endpoints_y[2*i-1] = initial_position[i, 2] - R[i]
		endpoints_y[2*i] = initial_position[i, 2] + R[i]
	end
	# Visualization settings
	circle_resolution = 201
	theta = LinRange(0, 2*pi, circle_resolution)
	function ball_plot!(x, y)
		return plot!(x, y, xlims = (0, L), ylims = (0, H), aspect_ratio =1, seriestype = [:shape], lw = 0, legend = false, framestyle = :none, grid = false, ticks = false, windowsize = (1200, 900), background_color_inside = :white, background_color_outside = :lightskyblue4, fillalpha = 0.6, linealpha = 0)
	end
	t = 0
	iteration_number = 0
	my_plot = plot()
	for i=1:N
		if i==1
			my_plot = plot()
		end
		ball_plot!(positions[i][1, 1] .+ R[i]*cos.(theta), positions[i][1, 2] .+ R[i]*sin.(theta))
	end
	display(my_plot)
	sleep(0.025)
	epsilon = 1e-7
	while(t+epsilon<T)
		if(t+dt>T)
			dt = T - t
		end
		iteration_number+= 1
		# update, then check forces for next iteration
		# update
		t += dt
		for i=1:N
			velocities[i][iteration_number+1, 1] = velocities[i][iteration_number, 1] + accelerations[i][iteration_number, 1]*dt
			velocities[i][iteration_number+1, 2] = velocities[i][iteration_number, 2] + accelerations[i][iteration_number, 2]*dt
			positions[i][iteration_number+1, 1] = positions[i][iteration_number, 1] + velocities[i][iteration_number, 1]*dt + accelerations[i][iteration_number, 1]*(dt^2)/2
			positions[i][iteration_number+1, 2] = positions[i][iteration_number, 2] + velocities[i][iteration_number, 2]*dt + accelerations[i][iteration_number, 2]*(dt^2)/2
			if i==1
				my_plot = plot()
			end
			endpoints_x[2*i-1] = positions[i][iteration_number+1, 1] - R[i]
			endpoints_x[2*i] = positions[i][iteration_number+1, 1] + R[i]
			endpoints_y[2*i-1] = positions[i][iteration_number+1, 2] - R[i]
			endpoints_y[2*i] = positions[i][iteration_number+1, 2] + R[i]
			ball_plot!(positions[i][iteration_number, 1] .+ R[i]*cos.(theta), positions[i][iteration_number, 2] .+ R[i]*sin.(theta))
		end
		display(my_plot)
		sleep(0.025)

		# check for collisions
		# check for wall-ball collisions
		for i=1:N
			r = R[i]
			m = M[i]
			x = positions[i][iteration_number+1, 1]
			y = positions[i][iteration_number+1, 2]
			v_x = velocities[i][iteration_number+1, 1]
			v_y = velocities[i][iteration_number+1, 2]
			#a_x = 0
			#a_y = 0
			# wall collisions : elastic
			if ((positions[i][iteration_number+1, 1]<R[i])&&(velocities[i][iteration_number+1, 1] <0) || ((L-positions[i][iteration_number+1, 1])<R[i])&&velocities[i][iteration_number+1, 1]>0) # side wall collision
				if verbose
					print("Wall-ball collision !\n")
				end
				accelerations[i][iteration_number + 1, 1] += -2*velocities[i][iteration_number+1, 1]/dt # perfectly elastic shock, impulse calculation
			end
			if ((positions[i][iteration_number+1, 2]<R[i])&&(velocities[i][iteration_number+1, 2] <0) || ((L-positions[i][iteration_number+1, 2])<R[i])&&velocities[i][iteration_number+1, 2]>0) # roof/ground collision
				if verbose
					print("Wall-ball collision !\n")
				end
				accelerations[i][iteration_number + 1, 2] += -2*velocities[i][iteration_number+1, 2]/dt # perfectly elastic shock, impulse calculation
			end
		end
		# check for ball-ball collision
		# collision relation is "anti-reflexive" and symmetric
		relative_p_x, PCS_x = sweep_n_prune(endpoints_x[p_x], indices_x[p_x])
		relative_p_y, PCS_y = sweep_n_prune(endpoints_y[p_y], indices_y[p_y])
		p_x = p_x[relative_p_x]
		p_y = p_y[relative_p_y]
		PCS = intersect(PCS_x, PCS_y)
		for pair in PCS
			i = pair[1]
			j = pair[2]
			overlapped = overlapping_check(positions[i][iteration_number+1, :], positions[j][iteration_number+1, :], R[i], R[j])
			getting_closer = getting_closer_check(positions[i][iteration_number+1, :], positions[j][iteration_number+1, :], velocities[i][iteration_number+1, :], velocities[j][iteration_number+1, :])
			if overlapped&&getting_closer
				if verbose
					print("Ball ", i, " and ball ", j, " just collided !\n")
				end
				# Hypothesis : forces only along normal to contact plane
				x_1 = positions[i][iteration_number+1, :]
				x_2 = positions[j][iteration_number+1, :]
				n = (x_2 - x_1)/norm(x_2 - x_1)
				v_1n = dot(velocities[i][iteration_number+1, :], n)
				v_2n = dot(velocities[j][iteration_number+1, :], n)
				a_1n = 2*M[j]/((M[i] + M[j])*dt)*(v_2n - v_1n)
				a_2n = 2*M[i]/((M[i] + M[j])*dt)*(v_1n - v_2n)
				accelerations[i][iteration_number + 1, :] += a_1n*n
				accelerations[j][iteration_number + 1, :] += a_2n*n
			end
		end

		#for i=1:(N-1)
			#for j=(i+1):N
				#overlapped = overlapping_check(positions[i][iteration_number+1, :], positions[j][iteration_number+1, :], R[i], R[j])
				#getting_closer = getting_closer_check(positions[i][iteration_number+1, :], positions[j][iteration_number+1, :], velocities[i][iteration_number+1, :], velocities[j][iteration_number+1, :])
				#if (overlapped && getting_closer)
					#print("Ball-ball collision !\n")
					## Hypothesis : forces only along normal to contact plane
					#x_1 = positions[i][iteration_number+1, :]
					#x_2 = positions[j][iteration_number+1, :]
					#n = (x_2 - x_1)/norm(x_2 - x_1)
					#v_1n = dot(velocities[i][iteration_number+1, :], n)
					#v_2n = dot(velocities[j][iteration_number+1, :], n)
					#a_1n = 2*M[j]/((M[i] + M[j])*dt)*(v_2n - v_1n)
					#a_2n = 2*M[i]/((M[i] + M[j])*dt)*(v_1n - v_2n)
					#accelerations[i][iteration_number + 1, :] .+= a_1n*n
					#accelerations[j][iteration_number + 1, :] .+= a_2n*n
				#end
			#end
		#end
	end
	if save_gif
		anim = @animate for iteration_number=1:number_iterations
			for i=1:N
				if i==1
					my_plot = plot()
				end
				ball_plot!(positions[i][iteration_number, 1] .+ R[i]*cos.(theta), positions[i][iteration_number, 2] .+ R[i]*sin.(theta))
			end
		end
		gif(anim, "anim_fps15.gif", fps = 15)
	end
end
