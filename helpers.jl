using LinearAlgebra
function generate_circles(N, L, H)
	" generate non-overlapping circles in box of width L and height H "
	# upper and lower bounds for the radius
	upper_bound = 1
	lower_bound = sqrt(L*H/(2*N))/4
	radii = lower_bound .+ (upper_bound-lower_bound)*rand(Float64, N)
	centers = rand(Float64, (N, 2))
	centers[:, 1] *= L
	centers[:, 2] *= H
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
	end
	return centers, radii
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
