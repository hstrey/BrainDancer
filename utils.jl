using Plots
using Statistics
using LinearAlgebra

genimg(p::AbstractArray) = p ./ first(findmax(p))

function getcoordinates(img)
	hcat(([i.I...] for i in findall(img .> 0))...)
end

function showsegments(segments)
	imgs = []
	for i in segments.segment_labels
		img = Gray.(map(j-> j == i ? segment_mean(segments,i) : 1, labels_map(segments)))
		pp = plot(img, title="S$i")	
		push!(imgs, pp)
	end
	plot(imgs...,layout=length(imgs),axis=nothing,framestyle=:none)
end

function fitellipse(xy::Matrix)
	# design matrix
	D = let x = xy[1,:], y = xy[2,:]
		[x.*x x.*y y.*y x y ones(size(xy,2))]
	end
	# scatter matrix
	S = D' * D
	# constraint matrix
	C = zeros(6,6)
	C[1, 3] = 2
	C[2, 2] = -1
	C[3, 1] = 2
	# solve eigensystem
	F = eigen(inv(S) * C)
	F.vectors[:, findmax(F.values) |> last]	
end

function fitellipsedirect(xy::Matrix)
	D1, D2 = let x = xy[1,:], y = xy[2,:]
		# quadratic part of the design matrix
		[x.*x x.*y y.*y],
		# linear part of the design matrix
		[x y ones(size(xy,2))]
	end
	# quadratic part of the scatter matrix
	S1 = D1' * D1
	# combined part of the scatter matrix
	S2 = D1' * D2
	# linear part of the scatter matrix
	S3 = D2' * D2

	T = -inv(S3) * S2' # for getting a2 from a1
	M = S1 + S2 * T # reduced scatter matrix
	M = [M[3, :]./2 -M[2, :] M[1, :]./2] # premultiply by inv(C1)
	
	# solve eigensystem
	F = eigen(M)
	evec = F.vectors
	# F.vectors[:, findmax(F.values) |> last]	
	cond = 4 * evec[1, :] .* evec[3, :] - evec[2, :].^2 # evaluate a’Ca
	a1 = evec[:, findall(cond .> 0)] # eigenvector for min. pos. eigenvalue	
	#T, M, cond, 
	[a1; T * a1] |> vec
end

function canonical(params::Vector)
	(A,B,C,D,E,F) = params
	c1 = B^2-4*A*C 
	c2 = 2(A*E^2 + C*D^2 - B*D*E + c1*F)
	c3 = sqrt((A-C)^2+B^2)
	a = -sqrt(c2*(A+C+c3))/c1
	b = -sqrt(c2*(A+C-c3))/c1
	x = (2*C*D-B*E)/(c1)
	y = (2*A*E-B*D)/(c1)
	θ = (C - A - c3)/B |> atan
	a, b, x, y, θ, rad2deg(θ)
end