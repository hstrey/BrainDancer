### A Pluto.jl notebook ###
# v0.19.9

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

# ╔═╡ 5b8d1c62-fc8f-11ec-202f-5f63a94bc6bf
begin
	import Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	Pkg.add("Colors")
	Pkg.add("Images")
	Pkg.add("Measures")
	Pkg.add("DataFrames")
	Pkg.add("Statistics")
	Pkg.add("RecursiveArrayTools")
	Pkg.develop(path=joinpath(@__DIR__, "..", "NIfTI.jl"))
	using Plots, CSV, NIfTI, Colors, Images, Measures
	using DataFrames
	using Statistics
	using RecursiveArrayTools
end

# ╔═╡ 6704f270-9d91-4d63-bf53-16e29c246605
phantom = niread("/Users/hstrey/Desktop/Phantom_talk/Phantom dataset/epi/epi.nii")

# ╔═╡ 4ff856d5-de11-4876-acd4-7e48d88fa2fb
# saving header for later use
phantom_head = phantom.header

# ╔═╡ 2c3bcfa6-1607-42cd-8ba3-53cf2f19c3fc
vsize = voxel_size(phantom.header)    # In mm

# ╔═╡ 6116b291-ad7b-4191-9ec6-ce2fa0938121
# function to display images in red
Red(x) = RGB(x,0,0)

# ╔═╡ 462eed79-97b3-441a-ac9a-c70616acb166
plot(Gray.(phantom[:,:,10,400] ./ findmax(phantom[:,:,10,400])[1]) ,aspect_ratio=1.0,axis = nothing,framestyle=:none,title="tenth slice",size=(400,450))

# ╔═╡ b52d0c52-164a-4c75-a138-923d02869023
md"""
Please provide the range for good slices
$(@bind s html"<input type=range min=1 max=28 value=1>")
$(@bind e html"<input type=range min=1 max=28 value=28>")
"""

# ╔═╡ 5d9b4de4-7207-4f63-a8a3-b85f7823d736
begin
	images = []
	for i in 1:28
		img = phantom[:,:,i,1]
		imgmax = findmax(img)[1]
		if (i>=s) && (i<=e)
			pp = plot(Gray.(img ./ imgmax) ,aspect_ratio=1.0,axis = nothing,framestyle=:none,title="$i")
		else
			pp = plot(Red.(img ./ imgmax) ,aspect_ratio=1.0,axis = nothing,framestyle=:none,title="$i")
		end
		push!(images,pp)
	end
	p = plot((images...),layout=length(images),margins=-1.8mm)
end

# ╔═╡ 03771ec3-008b-4185-85fe-c498a5e60e66
phantom_ok = phantom[2:end,:,s:e,:]

# ╔═╡ 1fce80c1-e747-4d19-997c-31400e8a864c
phantom_static = mean(phantom_ok, dims=4)[:,:,:,1]

# ╔═╡ c8c923e0-ab23-4d0c-b370-7f4018ae0b43
# using the original NIfTI header and changing the dimensions
phantom_head.dim = (3,83,84,12,1,1,1,1)

# ╔═╡ 4205ab34-d85c-4b61-8311-94e2b5663b5b
ni_static = NIVolume(phantom_head, phantom_static)

# ╔═╡ a00819fc-f642-49a3-b742-05d8fd4a9f70
niwrite("static.nii",ni_static)

# ╔═╡ fe493d0b-1ac3-4ea9-a44f-5af78784e47d
md"""
Here we estimate middle and radius by hand
"""

# ╔═╡ e11e83b1-4f70-4622-a6bf-519286ae18d5
md"""
Estimate center
$(@bind h html"<input type=range min=1 max=84 value=42>")horz
$(@bind v html"<input type=range min=1 max=84 value=42>")vert
$(@bind r html"<input type=range min=1 max=40 value=10>")radius
"""

# ╔═╡ b9389e5d-82f0-45e1-9697-4a9af900a0fd
(v,h,r)

# ╔═╡ abd4e302-5be2-43e1-90f6-44f1ae3198b6
md"""
lower threshold
$(@bind lt html"<input type=range min=0 max=20000 value=0>")
"""

# ╔═╡ 108254f4-12f8-478b-907f-43d996bfc6af
md"""
we can also try a more automated way by getting the center of mass and then
estimate the mask from there
"""

# ╔═╡ ef460365-6f17-43ba-9c3e-abdecb92ef0d
begin
	x_range = 1:84
	y_range = 1:83
	x_grid = x_range' .* ones(length(y_range))
	y_grid = ones(length(x_range))' .* y_range
end

# ╔═╡ 230f1a54-609e-4cf1-981c-174792780656
size(phantom_static)[3]

# ╔═╡ dadffda4-fa90-4468-b0e9-65591ee42cd6
begin
	mask_list = []
	for i in 1:size(phantom_static)[3]
		mask_tmp = phantom_static[:,:,i] .> lt
		cm_x_tmp = sum(x_grid .* phantom_static[:,:,i])/sum(phantom_static[:,:,i])
		cm_y_tmp = sum(y_grid .* phantom_static[:,:,i])/sum(phantom_static[:,:,i])
		for ix in 1:84
			for iy in 1:83
				if (ix-cm_x_tmp)^2+(iy-cm_y_tmp)^2 < 16^2
					mask_tmp[ix,iy] = 1
				end
			end
		end
		push!(mask_list, mask_tmp)
	end
end

# ╔═╡ 7f4bb963-5abd-4c44-ac6d-89bd6e4e5ca6
md"""
pick slice
$(@bind pick_slice html"<input type=range min=1 max=9 value=1>")slice
"""

# ╔═╡ 3b794916-b6dc-40e1-874d-ea3db6533cf4
begin
	plot(Gray.(phantom_static[:,:,pick_slice] ./ findmax(phantom_static[:,:,pick_slice])[1]),
		aspect_ratio=1.0,
		axis = nothing,
		framestyle=:none,
		title="first ok slice",
		size=(400,450))
	hline!([h],color=:red,label="horz $h")
	vline!([v],color=:green,label="vert $v")
	phi = 0:0.01:2π
	circle_x = r .* cos.(phi) .+ v
	circle_y = r .* sin.(phi) .+ h
	plot!(circle_x, circle_y, label="radius $r")
end

# ╔═╡ 9518856b-5eed-4a99-8e51-f2a08e4f8228
maxint = findmax(phantom_static[:,:,pick_slice])

# ╔═╡ b03fc810-65a8-43b9-a2a9-7051f1979e64
minint = findmin(phantom_static[:,:,pick_slice])

# ╔═╡ e7172e60-f796-45d0-9e63-c4aefe37a45a
mean(phantom_static[:,:,pick_slice])

# ╔═╡ c02d0ede-4b74-41ac-844c-8746532d14d2
mask = phantom_static[:,:,pick_slice] .> lt

# ╔═╡ 11e82358-e5e6-48d9-998a-7a01f1ebbb92
plot(Red.(mask))

# ╔═╡ 341eff6d-1a32-4aa3-a6f0-17bbe53997c4
cm_x = sum(x_grid .* phantom_static[:,:,pick_slice])/sum(phantom_static[:,:,pick_slice])

# ╔═╡ 5f1204b2-70be-42db-93a7-841a3bc81671
cm_y = sum(y_grid .* phantom_static[:,:,pick_slice])/sum(phantom_static[:,:,pick_slice])

# ╔═╡ bdd11d90-b2a4-4d3c-8647-e2e68a470cfc
plot(Red.(mask_list[pick_slice]))

# ╔═╡ 2f009527-88d5-4bd3-b3e6-ac756e8fb2c6
ml = Int16.(convert(Array,VectorOfArray(mask_list)))

# ╔═╡ 25e098c1-53fb-478c-93e2-5f4723610013
ni_mask = NIVolume(phantom_head, ml)

# ╔═╡ d25ff77c-3b2d-4f0a-a798-83f12643e142
niwrite("mask.nii",ni_mask)

# ╔═╡ Cell order:
# ╠═5b8d1c62-fc8f-11ec-202f-5f63a94bc6bf
# ╠═6704f270-9d91-4d63-bf53-16e29c246605
# ╠═4ff856d5-de11-4876-acd4-7e48d88fa2fb
# ╠═2c3bcfa6-1607-42cd-8ba3-53cf2f19c3fc
# ╠═6116b291-ad7b-4191-9ec6-ce2fa0938121
# ╠═462eed79-97b3-441a-ac9a-c70616acb166
# ╟─5d9b4de4-7207-4f63-a8a3-b85f7823d736
# ╠═b52d0c52-164a-4c75-a138-923d02869023
# ╠═03771ec3-008b-4185-85fe-c498a5e60e66
# ╠═1fce80c1-e747-4d19-997c-31400e8a864c
# ╠═c8c923e0-ab23-4d0c-b370-7f4018ae0b43
# ╠═4205ab34-d85c-4b61-8311-94e2b5663b5b
# ╠═a00819fc-f642-49a3-b742-05d8fd4a9f70
# ╟─fe493d0b-1ac3-4ea9-a44f-5af78784e47d
# ╠═3b794916-b6dc-40e1-874d-ea3db6533cf4
# ╠═e11e83b1-4f70-4622-a6bf-519286ae18d5
# ╠═b9389e5d-82f0-45e1-9697-4a9af900a0fd
# ╠═9518856b-5eed-4a99-8e51-f2a08e4f8228
# ╠═b03fc810-65a8-43b9-a2a9-7051f1979e64
# ╠═e7172e60-f796-45d0-9e63-c4aefe37a45a
# ╟─abd4e302-5be2-43e1-90f6-44f1ae3198b6
# ╠═11e82358-e5e6-48d9-998a-7a01f1ebbb92
# ╠═c02d0ede-4b74-41ac-844c-8746532d14d2
# ╠═108254f4-12f8-478b-907f-43d996bfc6af
# ╠═ef460365-6f17-43ba-9c3e-abdecb92ef0d
# ╠═230f1a54-609e-4cf1-981c-174792780656
# ╠═341eff6d-1a32-4aa3-a6f0-17bbe53997c4
# ╠═5f1204b2-70be-42db-93a7-841a3bc81671
# ╠═dadffda4-fa90-4468-b0e9-65591ee42cd6
# ╠═bdd11d90-b2a4-4d3c-8647-e2e68a470cfc
# ╠═7f4bb963-5abd-4c44-ac6d-89bd6e4e5ca6
# ╠═2f009527-88d5-4bd3-b3e6-ac756e8fb2c6
# ╠═25e098c1-53fb-478c-93e2-5f4723610013
# ╠═d25ff77c-3b2d-4f0a-a798-83f12643e142
