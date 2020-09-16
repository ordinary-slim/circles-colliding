using Plots
using LinearAlgebra
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
