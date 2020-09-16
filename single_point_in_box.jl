using Plots
using LinearAlgebra
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
