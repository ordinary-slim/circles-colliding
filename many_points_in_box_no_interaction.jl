using Plots
using LinearAlgebra

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
