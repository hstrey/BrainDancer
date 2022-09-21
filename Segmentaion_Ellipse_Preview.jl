### A Pluto.jl notebook ###
# v0.19.11

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
    import Pkg
    Pkg.activate(".")
	# for p in ["Plots", "Colors", "NIfTI", "Images", "ImageSegmentation", "ImageMorphology"]
	# 	Pkg.add(p)
	# end
	
	using Plots, Colors
	using NIfTI
	using Images, ImageSegmentation, ImageMorphology
	using Statistics
	using LinearAlgebra
end

# ╔═╡ d6ce949d-d0ef-4b9b-a7e1-b2a7c2258213
# Load auxiliary functionality
include("utils.jl")

# ╔═╡ 03fe33ed-4c31-4c3b-b402-494509eae5c2
const DATA_DIR = "../PhantomData"

# ╔═╡ 5a4c6be4-7657-47cd-ae84-dedd2f1facd9
# Load data
phantom, masks = niread(joinpath(DATA_DIR, "static.nii")),
 				 niread(joinpath(DATA_DIR, "mask.nii"));

# ╔═╡ 92d88bd4-c8e5-4167-bc40-011b44d3d021
# Setup visualization params
gr(display_type=:inline) 

# ╔═╡ b68e388a-2967-464f-9c90-cdb67eb97e2e
md"# Show images"

# ╔═╡ 919d2a18-1c0a-4d4d-a29c-7b50e10549fa
let images = []
	for i in 1:size(phantom,3)
		pp = plot(Gray.(genimg(phantom[:,:,i])),aspect_ratio=1.0,axis = nothing,framestyle=:none,title="$i")	
		push!(images,pp)
	end
	plot((images...),layout=length(images))
end

# ╔═╡ 1ca425e8-8a91-44f7-baa6-eed2c176e1ec
md"""
Preview parameters:
- Image #: $(@bind imgid html"<input type=range min=1 max=12 value=1></input>")
- Total segments: $(@bind segcount html"<input type=number value=4></input>")
- Segmentation threshold: $(@bind segcth html"<input type=number value=75></input>")
- Segment #: $(@bind segno html"<input type=number value=1></input>")
"""

# ╔═╡ 38e41da3-fdc1-4221-81c6-69de56162066
let img = genimg(phantom[:,:,imgid])
	# segment image
	segments = felzenszwalb(img, segcount, segcth)
	# select segment
	seg = Gray.(map(j-> j == segno ? 1 : 0, labels_map(segments)))
	# get segment border
	mseg = morpholaplace(seg)
	# estimate ellipse params
	E = mseg |> getcoordinates |> fitellipsedirect |> canonical
	# show center	
	img[round.(Int, E[3:4])...] = 1
	# show ellipse
	for (x,y) in [round.(Int, (E[3]-E[1]*sin(t), E[4]+E[2]*cos(t))) for t in 0:0.1:2π]
		img[x,y] = 1
	end
	# show image
	pseg = showsegments(segments)
	pimg = plot(Gray.(img) ,aspect_ratio=1.0,axis = nothing,framestyle=:none,title="img #$imgid",size=(300,350))
	plot(pimg, pseg)
end

# ╔═╡ Cell order:
# ╠═8d4361ce-39db-11ed-130c-7fb4d02c0f6e
# ╠═03fe33ed-4c31-4c3b-b402-494509eae5c2
# ╠═5a4c6be4-7657-47cd-ae84-dedd2f1facd9
# ╠═92d88bd4-c8e5-4167-bc40-011b44d3d021
# ╠═d6ce949d-d0ef-4b9b-a7e1-b2a7c2258213
# ╟─b68e388a-2967-464f-9c90-cdb67eb97e2e
# ╟─919d2a18-1c0a-4d4d-a29c-7b50e10549fa
# ╟─1ca425e8-8a91-44f7-baa6-eed2c176e1ec
# ╟─38e41da3-fdc1-4221-81c6-69de56162066
