### A Pluto.jl notebook ###
# v0.19.22

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
	Pkg.add("NIfTI")
	Pkg.add("Plots")
	Pkg.add("CSV")
	using Plots, CSV, NIfTI, Colors, Images, Measures
	using DataFrames
	using Statistics
	using RecursiveArrayTools
end

# ╔═╡ ef9548bd-042b-461d-a020-3dcee0efadf7
md"""
It looks like that there is always an empty row in the NIfTI file.  This row needs to be removed.  
"""

# ╔═╡ 6704f270-9d91-4d63-bf53-16e29c246605
phantom = niread("../PhantomData/UConn phantom data/104/104.nii")

# ╔═╡ 8f178d20-c6c2-4e55-8bf4-9cb016b785ee
size(phantom)

# ╔═╡ 4ff856d5-de11-4876-acd4-7e48d88fa2fb
# saving header for later use
phantom_head = phantom.header

# ╔═╡ 2c3bcfa6-1607-42cd-8ba3-53cf2f19c3fc
vsize = voxel_size(phantom.header)    # In mm

# ╔═╡ 6116b291-ad7b-4191-9ec6-ce2fa0938121
# function to display images in red
Red(x) = RGB(x,0,0)

# ╔═╡ 462eed79-97b3-441a-ac9a-c70616acb166
plot(Gray.(phantom[:,:,28,400] ./ findmax(phantom[:,:,28,400])[1]) ,aspect_ratio=1.0,axis = nothing,framestyle=:none,title="twentyeigth slice",size=(400,450))

# ╔═╡ 3b0a3f7c-2c64-46f6-a6ea-2577ba558a0c
md"""
The first step in the analysis is for the user to pick the good slices.  For this we want to display all slices and then the user can select the range of slices that are undistorted.  I see that we have to strike a compromise between displaying many images at the same time and their resolution
"""

# ╔═╡ b52d0c52-164a-4c75-a138-923d02869023
md"""
Please provide the range for good slices
$(@bind s html"<input type=range min=1 max=60 value=2>")
$(@bind e html"<input type=range min=1 max=60 value=59>")
"""

# ╔═╡ 5d9b4de4-7207-4f63-a8a3-b85f7823d736
begin
	images = []
	for i in 26:47
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
# for bias field correction we only need static ok images
# also disregard empty row
# stack_start
begin
	if s > 1
		stack_start = s-1
	else
		stack_start = 1
	end
	if e < size(phantom)[3]
		stack_end = e+1
	else
		stack_end = e
	end
	# only use first 200 time points - in first 200 time points phantom is not moving
	phantom_ok = phantom[2:end,:,stack_start:stack_end,1:200]
end

# ╔═╡ 43130725-c005-4d5f-b86e-c374f58311f4
aqui = CSV.read("../PhantomData/UConn phantom data/104/acq_times_104.csv",DataFrame)

# ╔═╡ 8549980d-f4d6-4be4-9ed5-ef4fc27481e8
p_log = CSV.read("../PhantomData/UConn phantom data/104/log104.csv",DataFrame)

# ╔═╡ a06ee4f8-31bb-4b7a-8c04-e4cf65ca53cc
max_motion = findmax(p_log[!,"Tmot"])[1]

# ╔═╡ 05ceb6d2-d3e0-4c3d-96bb-c38bb47dac0c
slices_without_motion = aqui[!,"Slice"][aqui[!,"Time"] .> max_motion]

# ╔═╡ 21b7d923-848d-4850-bfa2-7bd64102aad2
slices_ok = sort(slices_without_motion[s-1 .<= slices_without_motion .<= e+1])

# ╔═╡ 33472775-2fad-451a-8715-0bd5d4d87c53
slices_selected = collect((s-1):(e+1))

# ╔═╡ f21cb7c6-c0f9-49ce-a0b0-6ac95ee63523
slices_wm = [x in slices_ok ? 1 : 0 for x in slices_selected]

# ╔═╡ 3b43be1c-c782-4a82-9c76-fadce64b71e1
slices_df = DataFrame(Dict(:slice => slices_selected, :no_motion => slices_wm))

# ╔═╡ 200bf9f2-59d1-4503-90b3-8514198db6f6
CSV.write("../PhantomData/UConn phantom data/104/slices.csv",slices_df)

# ╔═╡ 1fce80c1-e747-4d19-997c-31400e8a864c
# average over the first 200 static slices
phantom_static = mean(phantom_ok, dims=4)[:,:,:,1]

# ╔═╡ c8c923e0-ab23-4d0c-b370-7f4018ae0b43
# using the original NIfTI header and changing the dimensions
phantom_head.dim = (3,83,84,stack_end-stack_start+1,1,1,1,1)

# ╔═╡ 94258efc-5f30-4d9d-8c1a-9710357866cf
ni_static = NIVolume(phantom_head, phantom_static)

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
	""
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
$(@bind pick_slice html"<input type=range min=1 max=47 value=1>")slice
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

# ╔═╡ b9389e5d-82f0-45e1-9697-4a9af900a0fd
(v,h,r,pick_slice)

# ╔═╡ 9518856b-5eed-4a99-8e51-f2a08e4f8228
maxint = findmax(phantom_static[:,:,pick_slice])

# ╔═╡ b03fc810-65a8-43b9-a2a9-7051f1979e64
minint = findmin(phantom_static[:,:,pick_slice])

# ╔═╡ e7172e60-f796-45d0-9e63-c4aefe37a45a
mean(phantom_static[:,:,pick_slice])

# ╔═╡ c02d0ede-4b74-41ac-844c-8746532d14d2
begin
	mask = phantom_static[:,:,pick_slice] .> lt
	"create mask"
end

# ╔═╡ 11e82358-e5e6-48d9-998a-7a01f1ebbb92
plot(Red.(mask))

# ╔═╡ 341eff6d-1a32-4aa3-a6f0-17bbe53997c4
cm_x = sum(x_grid .* phantom_static[:,:,pick_slice])/sum(phantom_static[:,:,pick_slice])

# ╔═╡ 5f1204b2-70be-42db-93a7-841a3bc81671
cm_y = sum(y_grid .* phantom_static[:,:,pick_slice])/sum(phantom_static[:,:,pick_slice])

# ╔═╡ bdd11d90-b2a4-4d3c-8647-e2e68a470cfc
plot(Red.(mask_list[pick_slice]))

# ╔═╡ 2f009527-88d5-4bd3-b3e6-ac756e8fb2c6
ml = Float32.(convert(Array,VectorOfArray(mask_list)))

# ╔═╡ 25e098c1-53fb-478c-93e2-5f4723610013
ni_mask = NIVolume(phantom_head, ml)

# ╔═╡ 9fcb4683-d76c-4f15-b3a5-63d0b2ba9dcc
Good_slices = phantom[2:end,:,stack_start:stack_end,:]

# ╔═╡ f44507ee-501d-4a7d-8019-b314272b90bd
phantom_slices = phantom.header

# ╔═╡ 613decfa-c2c3-4426-af1a-4fb68e5839b5
phantom_slices.dim = (4,83,84,stack_end-stack_start+1,800,1,1,1)

# ╔═╡ d5cecab3-1d59-43bf-92fb-4084247a8e85
stack_end-stack_start+1

# ╔═╡ 3afa8868-6834-4fd2-a76a-065ea531b3f1
ni_whole = NIVolume(phantom_slices, Good_slices)

# ╔═╡ d25ff77c-3b2d-4f0a-a798-83f12643e142
# write out the static phantom and its mask
begin
	niwrite("../PhantomData/UConn phantom data/104/mask.nii",ni_mask)
	niwrite("../PhantomData/UConn phantom data/104/static.nii",ni_static)
	niwrite("../PhantomData/UConn phantom data/104/Good_slices.nii",ni_whole)
end

# ╔═╡ 87b6f9b0-7a4f-4e0c-99ce-1f5adb65ce60
# no verbose; if want verbose, -v
ANTscommand = `N4BiasFieldCorrection -s 1 -d 3 -b \[ 200,3 \] -i "../PhantomData/UConn phantom data/104/static.nii" -x ../PhantomData/UConn phantom data/104/mask.nii -o \[ ../PhantomData/UConn phantom data/104/BFC_static200.nii, ../PhantomData/UConn phantom data/104/BFC_bias200.nii \]`

# ╔═╡ 296d707f-a65c-4e21-becb-f6dc06e2b3e8
typeof(ANTscommand)

# ╔═╡ 852d086f-9b13-4f7c-b831-1dd2d68378b9
# ╠═╡ disabled = true
#=╠═╡
run(ANTscommand)
  ╠═╡ =#

# ╔═╡ 6acaf6cc-c43c-4ba6-97f1-7842129a3d38
bias = niread("../PhantomData/UConn phantom data/104/BFC_bias200.nii")

# ╔═╡ 7e895010-048d-4fde-ad9c-684b979319a4
Corrected = niread("../PhantomData/UConn phantom data/104/BFC_static200.nii")

# ╔═╡ f66b5dd2-5b16-425d-830d-7006bdac3bde
begin
	images2 = []
	for i in 1:size(Corrected)[3]
		img = Corrected[:,:,i]
		imgmax = findmax(img)[1]
		pp = plot(Gray.(img ./ imgmax) ,aspect_ratio=1.0,axis = nothing,framestyle=:none,title="$i")
		push!(images2,pp)
	end
	p2 = plot((images2...),layout=length(images2),margins=-1.8mm)
end

# ╔═╡ b6e64653-42ef-492b-8489-466fd4f3c6f1
begin 
	BFC_image = zeros(size(ni_whole))
	for i in 1:size(ni_whole)[4]
		for j in 1:size(ni_whole)[3]
			BFC_image[:,:,j,i] = ni_whole[:,:,j,i]./bias[:,:,j]
		end
	end
end

# ╔═╡ 818884d8-2e0a-468d-86b3-900790bf0d4d
md"""
Pick a time slice
$(@bind t html"<input type=number min=1 max=800 value=42>")
"""

# ╔═╡ 5b21716f-0cbf-48cb-aacc-a21c49dc81bf
begin
	images1 = []
	for i in 1:e-s+1
		img1 = BFC_image[:,:,i,t]
		imgmax1 = findmax(img1)[1]
		
		pp1 = plot(Gray.(img1 ./ imgmax1) ,aspect_ratio=1.0,axis = nothing,framestyle=:none,title="$i")
		
		push!(images1,pp1)
	end
	p1 = plot((images1...),layout=length(images1),margins=-1.8mm)
end

# ╔═╡ b5c2f9bb-f487-4a92-904d-a8714eb0da04
BFC_NIfTI = NIVolume(phantom_slices, BFC_image)

# ╔═╡ 501049c5-239b-4b16-bd1b-406d5ed96ab2
niwrite("../PhantomData/UConn phantom data/104/BFC_time_series.nii",BFC_NIfTI)

# ╔═╡ Cell order:
# ╠═5b8d1c62-fc8f-11ec-202f-5f63a94bc6bf
# ╠═ef9548bd-042b-461d-a020-3dcee0efadf7
# ╠═6704f270-9d91-4d63-bf53-16e29c246605
# ╠═8f178d20-c6c2-4e55-8bf4-9cb016b785ee
# ╠═4ff856d5-de11-4876-acd4-7e48d88fa2fb
# ╠═2c3bcfa6-1607-42cd-8ba3-53cf2f19c3fc
# ╠═6116b291-ad7b-4191-9ec6-ce2fa0938121
# ╠═462eed79-97b3-441a-ac9a-c70616acb166
# ╠═3b0a3f7c-2c64-46f6-a6ea-2577ba558a0c
# ╟─5d9b4de4-7207-4f63-a8a3-b85f7823d736
# ╠═b52d0c52-164a-4c75-a138-923d02869023
# ╠═03771ec3-008b-4185-85fe-c498a5e60e66
# ╠═43130725-c005-4d5f-b86e-c374f58311f4
# ╠═8549980d-f4d6-4be4-9ed5-ef4fc27481e8
# ╠═a06ee4f8-31bb-4b7a-8c04-e4cf65ca53cc
# ╠═05ceb6d2-d3e0-4c3d-96bb-c38bb47dac0c
# ╠═21b7d923-848d-4850-bfa2-7bd64102aad2
# ╠═33472775-2fad-451a-8715-0bd5d4d87c53
# ╠═f21cb7c6-c0f9-49ce-a0b0-6ac95ee63523
# ╠═3b43be1c-c782-4a82-9c76-fadce64b71e1
# ╠═200bf9f2-59d1-4503-90b3-8514198db6f6
# ╠═1fce80c1-e747-4d19-997c-31400e8a864c
# ╠═c8c923e0-ab23-4d0c-b370-7f4018ae0b43
# ╠═94258efc-5f30-4d9d-8c1a-9710357866cf
# ╟─fe493d0b-1ac3-4ea9-a44f-5af78784e47d
# ╠═3b794916-b6dc-40e1-874d-ea3db6533cf4
# ╠═e11e83b1-4f70-4622-a6bf-519286ae18d5
# ╠═b9389e5d-82f0-45e1-9697-4a9af900a0fd
# ╠═9518856b-5eed-4a99-8e51-f2a08e4f8228
# ╠═b03fc810-65a8-43b9-a2a9-7051f1979e64
# ╠═e7172e60-f796-45d0-9e63-c4aefe37a45a
# ╟─abd4e302-5be2-43e1-90f6-44f1ae3198b6
# ╠═11e82358-e5e6-48d9-998a-7a01f1ebbb92
# ╟─c02d0ede-4b74-41ac-844c-8746532d14d2
# ╠═108254f4-12f8-478b-907f-43d996bfc6af
# ╟─ef460365-6f17-43ba-9c3e-abdecb92ef0d
# ╠═230f1a54-609e-4cf1-981c-174792780656
# ╠═341eff6d-1a32-4aa3-a6f0-17bbe53997c4
# ╠═5f1204b2-70be-42db-93a7-841a3bc81671
# ╠═dadffda4-fa90-4468-b0e9-65591ee42cd6
# ╠═bdd11d90-b2a4-4d3c-8647-e2e68a470cfc
# ╠═7f4bb963-5abd-4c44-ac6d-89bd6e4e5ca6
# ╠═2f009527-88d5-4bd3-b3e6-ac756e8fb2c6
# ╠═25e098c1-53fb-478c-93e2-5f4723610013
# ╠═9fcb4683-d76c-4f15-b3a5-63d0b2ba9dcc
# ╠═f44507ee-501d-4a7d-8019-b314272b90bd
# ╠═613decfa-c2c3-4426-af1a-4fb68e5839b5
# ╠═d5cecab3-1d59-43bf-92fb-4084247a8e85
# ╠═3afa8868-6834-4fd2-a76a-065ea531b3f1
# ╠═d25ff77c-3b2d-4f0a-a798-83f12643e142
# ╠═87b6f9b0-7a4f-4e0c-99ce-1f5adb65ce60
# ╠═296d707f-a65c-4e21-becb-f6dc06e2b3e8
# ╠═852d086f-9b13-4f7c-b831-1dd2d68378b9
# ╠═6acaf6cc-c43c-4ba6-97f1-7842129a3d38
# ╠═7e895010-048d-4fde-ad9c-684b979319a4
# ╠═f66b5dd2-5b16-425d-830d-7006bdac3bde
# ╠═b6e64653-42ef-492b-8489-466fd4f3c6f1
# ╠═5b21716f-0cbf-48cb-aacc-a21c49dc81bf
# ╠═818884d8-2e0a-468d-86b3-900790bf0d4d
# ╠═b5c2f9bb-f487-4a92-904d-a8714eb0da04
# ╠═501049c5-239b-4b16-bd1b-406d5ed96ab2
