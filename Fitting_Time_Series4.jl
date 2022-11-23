### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8d4361ce-39db-11ed-130c-7fb4d02c0f6e
begin
    # import Pkg
    # Pkg.activate(".")
	# for p in ["Plots", "Colors", "NIfTI", "Images", "ImageSegmentation", "ImageTransformations", "LsqFit", "CoordinateTransformations", "Rotations"]
	# 	Pkg.add(p)
	# end
	
	using Plots, Colors, StatsPlots, Measures
	using NIfTI
	using Images, ImageSegmentation, ImageTransformations
	using CoordinateTransformations, Rotations
	using Statistics
	using LinearAlgebra, LsqFit, Interpolations
	using DataFrames, CSV
	const ITI = ImageTransformations.Interpolations
end

# ╔═╡ d6ce949d-d0ef-4b9b-a7e1-b2a7c2258213
# Load auxiliary functionality
include("utils.jl")

# ╔═╡ 03fe33ed-4c31-4c3b-b402-494509eae5c2
const DATA_DIR = "PhantomData"

# ╔═╡ f58ba902-0e84-435c-ab41-1ce0c5660830
# Load data
phantom_ts = niread(joinpath(DATA_DIR, "BFC_time_series.nii"));

# ╔═╡ 15ed4afa-3339-4938-b647-a6403ab90078
sz = size(phantom_ts)

# ╔═╡ ac243c5f-dcc8-4b1e-b63f-0cb96bb6f6d6
# valid range for the volume
zs = 1:sz[3] #-1

# ╔═╡ 92d88bd4-c8e5-4167-bc40-011b44d3d021
# Setup visualization params
gr(display_type=:inline) 

# ╔═╡ b68e388a-2967-464f-9c90-cdb67eb97e2e
md"**Show averaged slices**"

# ╔═╡ 82544c78-379f-4bfa-8f70-91a92b40fdce
# volume averaged over 200 static slices
staticimgs = dropdims(mean(view(phantom_ts,:,:,:,1:200), dims=4), dims=4);

# ╔═╡ 919d2a18-1c0a-4d4d-a29c-7b50e10549fa
let images = []
	for i in 1:size(staticimgs,3)
		pp = plot(Gray.(genimg(staticimgs[:,:,i])),aspect_ratio=1.0,axis = nothing,framestyle=:none,title="$i")	
		push!(images,pp)
	end
	plot((images...),layout=length(images))
end

# ╔═╡ 006cc58c-411c-4fb2-b219-45a86df742d8
md"**Fit ellipse for *every* averaged Z-slice**"

# ╔═╡ 868be102-6d42-4e71-94df-c8486558cc53
# get ellipse parameters over averaged statis volume
staticEs = let res = zeros(9,length(zs))	
	for i in zs
		res[:,i] = getellipse(view(staticimgs,:,:,i), verbose=false)
	end
	res
end

# ╔═╡ c3b152da-97da-42b3-93e7-a769a9b751c3
plot(staticEs[1,:], label="x0")

# ╔═╡ 6fed789f-c68f-48d7-8ef8-a10d1dfce7dd
plot(staticEs[2,:], label="y0")

# ╔═╡ 497edfd4-8c69-41cd-a9fb-35a3cf35c23b
md"**Fit line through the slices' centers & determine a deviation angle from `z`-axis**"

# ╔═╡ 4e003878-246b-46ff-8496-eff08b3a0596
lineparams = let data = staticEs[1:2,2:end-1] # use slices from 2 to 12
	μ = mean(data, dims=2)
	F = svd(data .- μ)
	dir = vec(F.U[1,:])
	vec(μ), dir, dir[2]/dir[1] # mean, direction, slope
end

# ╔═╡ 60657161-9d68-43b8-8f69-20489a76f786
md"**Calculate transformation parameters for phantom data generation**"

# ╔═╡ 4a1290c6-ea27-4851-b4e6-3a8956b527b4
# A deviation angle from z-axis
θ = let r = π/2-atan(lineparams[3])
	@info "Degree" d=rad2deg(r)
	r
end

# ╔═╡ fc457183-7e42-4da7-96c0-34815234a8da
# A rotation angle within xy-plain
αs = staticEs[5,:]

# ╔═╡ 2e28a4f2-096a-4257-87cf-e91ef9469786
# An adjusted ellipse centers of the each slice
centers = let α = mean(αs),
			  rng = -1:0.15:1.0,
		      cc = staticEs[1:2,:],
			  # xc = (x0,z)->z.*cos(θ).*cos(αs[z]).+x0,
			  # yc = (y0,z)->z.*cos(θ).*sin(αs[z]).+y0
		      xc = (x0,z)->z.*cos(θ).*cos(α).+x0,
			  yc = (y0,z)->z.*cos(θ).*sin(α).+y0
	# hcat([[xc(x,z/14), yc(y,z/14)] for (z,(x,y)) in enumerate(eachcol(cc))]...)
	hcat([[xc(x,z), yc(y,z)] for (z,(x,y)) in zip(rng, eachcol(cc))]...)	
	# xy = lineparams[1] .+ lineparams[2].*rng'
end

# ╔═╡ 5416e834-2d41-4c5c-b0bf-dadb9eb9af3c
let x = staticEs[1,:], y=staticEs[2,:], rng=collect(-1:0.15:1.)
	p = scatter(x, y, label="centers", legend=:bottomright)
	xy = lineparams[1] .+ lineparams[2].*rng'
	plot!(p, xy[1,:], xy[2,:], label="fit")	
	plot!(p, centers[1,:], centers[2,:], label="sim")
end

# ╔═╡ e6deda84-ff89-4998-852d-cbf752d98936
md"**Define rotation function for coordinate transformation given the ellipse parameters from the static slices**"

# ╔═╡ 66e9f543-3c0a-4344-a067-3fbef6228f1d
"""
	rotate(x, y, z, a, b, c, u, v, w, θ)

Rotate point `(x,y,z)` around the line that passes through the point `(a,b,c)` in the direction of `<u,v,w>` on the angle `θ`.
"""
function rotate(x, y, z, a, b, c, u, v, w, θ)
	(a*(v^2+w^2)−u*(b*v+c*w−u*x−v*y−w*z))*(1−cos(θ))+x*cos(θ)+(−c*v+b*w−w*y+v*z)*sin(θ),
	(b*(u^2+w^2)−v*(a*u+c*w−u*x−v*y−w*z))*(1−cos(θ))+y*cos(θ)+(c*u−a*w+w*x−u*z)*sin(θ),
	(c*(u^2+v^2)−w*(a*u+b*v−u*x−v*y−w*z))*(1−cos(θ))+z*cos(θ)+(−b*u+a*v−v*x+u*y)*sin(θ)
end

# ╔═╡ 8bb2090f-6ef8-445b-bca8-e124293bc459
"""
	rotatevoxel(o::Vector, p::Vector) -> f(x,y,z,θ)

Return a function that perform rotation of the voxel with coordinates `(x,y,z)` on the angle `θ` within XY-plane around the line started in an origin point `o` and passing through a point `p`.
"""
function rotatevoxel(origin, p)
	dir = p .- origin |> normalize	
	(x,y,z,θ) -> rotate(x, y, z, origin..., dir..., θ)
end

# ╔═╡ 3ff8e722-ec91-4751-bff6-e076599d2855
cc = let iidxs = (43, 43, 1), θ = staticEs[5],
		 origin = [staticEs[:,1][1:2]; 1],
		 last = [staticEs[:,2][1:2]; sz[3]]
	rvfn = rotatevoxel(origin, last) # create rotation function
	cc = rvfn(iidxs..., θ)       # rotate coordinates
	iidxs, θ, cc
end

# ╔═╡ 7bba97a1-9850-4dd1-a740-37edcdf44e1b
md"**Create an interpolation function from the static phantom avarage**"

# ╔═╡ 81b515d0-f00c-457a-babc-2480692ad8e7
phantom_itp = let (r,c,h) = size(staticimgs),
				  xs = 1:r, ys = 1:c, zs = 1:h
	extrapolate(
		scale(interpolate(staticimgs, BSpline(Linear())), xs, ys, zs),
		Line()
	) 
end;

# ╔═╡ cf938c35-3843-4fe2-a6f8-be3381f39141
# Test interpolation
let actcc = cc[3], intcc = round.(Int, actcc)
	v1 = staticimgs[intcc...]
	@info "Actual" intcc v1
	v2 = phantom_itp(actcc...)
	@info "Interpolated" actcc v2
end

# ╔═╡ 0a597747-f5cd-40c5-a9a4-7005a34f6b94
md"**Load rotation data of dynamic phase**"

# ╔═╡ 2b3701ed-c053-4658-b693-f58cd565d699
df = CSV.read(joinpath(DATA_DIR, "epi",  "log.csv"), DataFrame);

# ╔═╡ c33d2d8f-abc9-4186-9b33-0af1a3508455
firstrotidx = findfirst(e->e>20, df.CurPos)

# ╔═╡ b1e37158-a0ca-4556-ab11-4beafd0e0ee5
# rotation angles in radians
angles = let quant = 2^13	
	[a > π ? a-2π : a  for a in (df.CurPos[firstrotidx:end] ./ quant).*(2π)]
end

# ╔═╡ 4dab67c1-b204-42fb-a621-06db4fb4ae11
plot(angles, ylab="α")

# ╔═╡ 1557be9b-2a83-45a1-9b81-ce72f27ba64c
md"**Generate intensities from the static avarage by rotating to an angle of a rotation of a dynamic slice**"

# ╔═╡ d59ccb58-f737-4b46-be43-ab6cc46679cc
md"""
Coordinates:
- X: $(@bind x html"<input type=number min=1 max=83 value=42></input>")
- Y: $(@bind y html"<input type=number min=1 max=84 value=42></input>")
- Z: $(@bind z html"<input type=number min=1 max=13 value=1></input>")
"""

# ╔═╡ 39ddb2af-3643-45ff-b8f0-f8c8e09bcbc7
let θs = angles, # rotation angles of dynamic volumes
	rvfn = rotatevoxel([centers[:,2]; 1], [centers[:,12]; 13]) # create rotation function		
	sim = [phantom_itp(rvfn(x, y, z, θ)...) for (i,θ) in enumerate(θs)]
	p = plot(sim, label="sim($x, $y, $z)", legend=:topright)
	acc = phantom_ts[x,y,z,firstrotidx:end]
	plot!(p, acc, label="actual")
	mean_sim = mean(sim)
	mean_acc = mean(acc)
	plot!(p, [1,length(θs)],[mean_sim,mean_sim], label="mean sim")
	plot!(p, [1,length(θs)],[mean_acc,mean_acc], label="mean acc")
end

# ╔═╡ ea1bcfbb-a214-4acc-9242-08fd1d867d75
begin
	@info "Actual" intencity=staticimgs[x,y,z]
	@info "Interpolated" intencity=phantom_itp(x,y,z)
end

# ╔═╡ 5901e7fc-9c1d-46f1-8c93-128106e974a5
md"**Generate simulated image**"

# ╔═╡ a56278b4-9b95-4a5a-aaa8-3e15a83db5bc
"""
	simulated_coordinates(sz::Tuple, a::Vector, b::Vector, θ::Float64)

Generate simulated coordinate set of size `sz` by rotating voxels within
*xy*-plain on angle `θ` around the line passing through points `a` and `b`.
"""
function simulated_coordinates(sz::Tuple, a, b, θ::Float64)
	rvfn = rotatevoxel(a, b) # create rotation function
	[rvfn(i,j,k,θ) for i in 1:sz[1], j in 1:sz[2], k in 1:sz[3]]
end

# ╔═╡ 6b16b27d-0dfe-4900-bcf5-efa806ce509c
"""
	simulated_coordinates_at_z(sz::Tuple, z, a::Vector, b::Vector, θ::Float64)

Generate simulated coordinate set of size `sz` at depth `z`, by rotating voxels within
*xy*-plain on angle `θ` around the line passing through points `a` and `b`.
"""
function simulated_coordinates_at_z(sz::Tuple, z, a, b, θ::Float64)
	rvfn = rotatevoxel(a, b) # create rotation function
	[rvfn(i,j,z,θ) for i in 1:sz[1], j in 1:sz[2]]
end

# ╔═╡ bacd2be5-c44b-4c2e-9660-7d7890ac7a01
md"""
Generate image from an avarage volume slice `Z` by rotating it at an angle `θ` ∈ [-π,π].

Slice (Z): $(@bind sliceId html"<input type=number min=1 max=13 value=1></input>")
Rotation angle (θ): $(@bind theta html"<input type=number min=-180 max=180 value=0></input>")
"""

# ╔═╡ 11841d64-1d44-4886-9b0e-7eadaa2b8da8
let θ = deg2rad(theta)
	# cc = staticEs[1:2,:]          # ellipese centers from static average
	cc = centers                    # simulated centers
	a = [cc[:,2]; 1]                # line from static average ellipse centers
	b = [cc[:,end-1]; 13]
	# generate image
	coords = simulated_coordinates_at_z(sz, sliceId, a, b, θ)
	sim = map(c->phantom_itp(c...), coords)
	gen = Gray.(sim |> genimg)
	# show averaged image	
	ave = Gray.(staticimgs[:,:,sliceId] |> genimg)
	pave = plot(ave, aspect_ratio=1.0, axis=nothing, framestyle=:none, title="img z=$sliceId", size=(300,350))
	# show generated image
	pgen = plot(gen, aspect_ratio=1.0, axis=nothing, framestyle=:none, title="generated at $theta")
	plot(pave, pgen)	
end

# ╔═╡ Cell order:
# ╠═8d4361ce-39db-11ed-130c-7fb4d02c0f6e
# ╠═03fe33ed-4c31-4c3b-b402-494509eae5c2
# ╠═f58ba902-0e84-435c-ab41-1ce0c5660830
# ╠═15ed4afa-3339-4938-b647-a6403ab90078
# ╠═ac243c5f-dcc8-4b1e-b63f-0cb96bb6f6d6
# ╠═92d88bd4-c8e5-4167-bc40-011b44d3d021
# ╠═d6ce949d-d0ef-4b9b-a7e1-b2a7c2258213
# ╟─b68e388a-2967-464f-9c90-cdb67eb97e2e
# ╠═82544c78-379f-4bfa-8f70-91a92b40fdce
# ╟─919d2a18-1c0a-4d4d-a29c-7b50e10549fa
# ╟─006cc58c-411c-4fb2-b219-45a86df742d8
# ╠═868be102-6d42-4e71-94df-c8486558cc53
# ╠═c3b152da-97da-42b3-93e7-a769a9b751c3
# ╠═6fed789f-c68f-48d7-8ef8-a10d1dfce7dd
# ╟─497edfd4-8c69-41cd-a9fb-35a3cf35c23b
# ╠═4e003878-246b-46ff-8496-eff08b3a0596
# ╟─60657161-9d68-43b8-8f69-20489a76f786
# ╟─4a1290c6-ea27-4851-b4e6-3a8956b527b4
# ╠═fc457183-7e42-4da7-96c0-34815234a8da
# ╠═2e28a4f2-096a-4257-87cf-e91ef9469786
# ╠═5416e834-2d41-4c5c-b0bf-dadb9eb9af3c
# ╟─e6deda84-ff89-4998-852d-cbf752d98936
# ╟─66e9f543-3c0a-4344-a067-3fbef6228f1d
# ╟─8bb2090f-6ef8-445b-bca8-e124293bc459
# ╠═3ff8e722-ec91-4751-bff6-e076599d2855
# ╟─7bba97a1-9850-4dd1-a740-37edcdf44e1b
# ╠═81b515d0-f00c-457a-babc-2480692ad8e7
# ╠═cf938c35-3843-4fe2-a6f8-be3381f39141
# ╟─0a597747-f5cd-40c5-a9a4-7005a34f6b94
# ╠═2b3701ed-c053-4658-b693-f58cd565d699
# ╠═c33d2d8f-abc9-4186-9b33-0af1a3508455
# ╠═b1e37158-a0ca-4556-ab11-4beafd0e0ee5
# ╟─4dab67c1-b204-42fb-a621-06db4fb4ae11
# ╟─1557be9b-2a83-45a1-9b81-ce72f27ba64c
# ╟─d59ccb58-f737-4b46-be43-ab6cc46679cc
# ╠═39ddb2af-3643-45ff-b8f0-f8c8e09bcbc7
# ╠═ea1bcfbb-a214-4acc-9242-08fd1d867d75
# ╟─5901e7fc-9c1d-46f1-8c93-128106e974a5
# ╟─a56278b4-9b95-4a5a-aaa8-3e15a83db5bc
# ╟─6b16b27d-0dfe-4900-bcf5-efa806ce509c
# ╟─bacd2be5-c44b-4c2e-9660-7d7890ac7a01
# ╠═11841d64-1d44-4886-9b0e-7eadaa2b8da8
