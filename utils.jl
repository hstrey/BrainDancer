using Statistics
using LinearAlgebra
using LsqFit
using Images
using ImageSegmentation

maxvalue(d) = reduce((x, y) -> d[x] = d[y] ? x : y, keys(d))

genimg(p::AbstractArray) = p ./ maximum(p)

function genunitimg!(img, arr)
    minval, maxval = extrema(arr)
    diff = maxval-minval
    broadcast!(x-> (x-minval)/diff, img, arr)
end

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
    x, y, a, b, θ, rad2deg(θ)
end

function gaussellipse(xy,p)
    x0,y0,rx,ry,θ,a,A,bg,σ = p #unpack parameters
    x = xy[:,1]
    y = xy[:,2]
    dx = x .- x0
    dy = y .- y0
    ct = cos(θ)
    st = sin(θ)
    return bg .- A * exp.( -(1 .-sqrt.(( dx .* ct .+ dy .* st ).^2/rx^2+(dx .* st .- dy .* ct).^2/ry^2)).^2/σ^2)
end

function fitellipselsq(xy::AbstractMatrix, z::Vector, p0)
    fit = curve_fit(gaussellipse, xy, z, p0)
    fit.param
end

function segment3(img::AbstractMatrix{T}) where {T}
    # outer boundary segmentation of the image
    sz = size(img)
    out_seeds = [(CartesianIndex(1,1), 1),
                 (CartesianIndex(sz[1],sz[2]), 1),
                 (CartesianIndex(sz[1]>>1,sz[2]>>1), 2)]
    outsegs = seeded_region_growing(img, out_seeds)

    # inneter boundary segmentation of the image
    idxs = findall(labels_map(outsegs) .== 2)
    in_seeds = [(idxs[1], 1), (idxs[10], 1), (idxs[end-10], 1), (idxs[end], 1),
                (CartesianIndex(sz[1]>>1,sz[2]>>1), 2)]
    insegs = seeded_region_growing(img, in_seeds)

    # combine segments
    push!(outsegs.segment_labels, 3)
    outsegs.segment_means[3] = segment_mean(insegs, 2)
    outsegs.segment_pixel_count[2] -= insegs.segment_pixel_count[2]
    outsegs.segment_pixel_count[3] = insegs.segment_pixel_count[2]
    outsegs.image_indexmap[insegs.image_indexmap .== 2] .= 3
    return outsegs
end

function edge3(segs::SegmentedImage)
    # get inner region
    inner = labels_map(segs) .== 3
    # get an inner edge
    img_edges = canny(inner, (Percentile(80), Percentile(20)), 2);
    # remove any points outside of outer boundary
    for ci in findall(labels_map(segs) .== 1)
        img_edges[ci] = 0
    end
    img_edges
end

"""
    prepare_initial_point(p::Vector)

Create an initial solution for fitting optimization problem
"""
function prepare_initial_point(p::Vector)
    x0, y0, a, b, θ, α, A, bg, σ = p
    # flip axes so a > b
    dont_flip = a > b
    a, b = dont_flip ? (a,b) : (b,a)
    # if the axes were fliiped rotate ellipse by π/2
    θ = dont_flip ? θ : (θ - π/2)
    #  make sure that angle within [0;π]
    θ %= π
    θ = θ < 0 ? π + θ : θ
    # form initial solution
    [x0, y0, a, b, θ, α, (A < 0 ? 0.0 : A), bg, abs(σ)]
end

function fitellipse3(img::AbstractMatrix{T}, segs::SegmentedImage, edge::Matrix{Bool};
                     verbose=true, keepinitialonerror=true, secondfit=true) where {T}
    r,c = size(img)
    # find ecllipse in edge
    coords = hcat(([i.I...] for i in findall(edge))...)
    E1 = fitellipsedirect(coords) |> canonical
    verbose && @debug "Initial fit" E1

    # get all points within the outer boundry
    coords2 = hcat(([ci.I...] for ci in findall(labels_map(segs) .!= 1))...)
    z = [img[x,y] for (x,y) in eachcol(coords2)]

    # refine ellipse paramaters using outer segment points
    # and edge ellipse estimate
    p0 = prepare_initial_point([E1[1:5]..., 0, extrema(z)..., 0.01])
    fit = curve_fit(gaussellipse, coords2', z, p0, autodiff=:forwarddiff)
    E = fit.param
    verbose && @debug "LSQ fit 1" p0 E

    # refine again
    lb = [c/2-3.0, r/2-3.0, 0.0, 0.0, -π, -π/2,   0,   0, 0.001]
    ub = [c/2+3.0, r/2+3.0, Inf, Inf,  π,  π/2, Inf, Inf, 1.000]
    p0 = prepare_initial_point(E)
    E = try
        if secondfit
            fit = curve_fit(gaussellipse, coords2', z, p0, lower=lb, upper=ub, autodiff=:forwarddiff)
            fit.param
        else
            p0
        end
    catch ex
        p0 = prepare_initial_point([E1[1:5]..., 0, minimum(z), maximum(z)/2, 0.01])
        if keepinitialonerror
            p0
        else
            @debug "Error. Trying with constraints." p0
            fit = curve_fit(gaussellipse, coords2', z, p0, lower=lb, upper=ub)
            fit.param
        end
    end
    verbose && secondfit && @debug "LSQ fit 2" p0 E

    return E
end

function getellipse(img::AbstractArray; verbose=true)
    segs = segment3(img)
    fitellipse3(img, segs, edge3(segs); verbose)
end

function gaussellipse3d(xyz, p)
    x0,y0,rx,ry,θ,α,A,bg,σ = p #unpack parameters
	ct, st = cos(θ), sin(θ)
	c1, c2 = ct*cos(α), ct*sin(α)
	# xc = z->z.*cos(θ).*cos(α).+x0
	# yc = z->z.*cos(θ).*sin(α).+y0
	# dx = map(r->r[1]-r[3]*c1+x0, eachrow(xyz))
	# dy = map(r->r[2]-r[3]*c2+y0, eachrow(xyz))
	# dx = view(xyz, :, 1) .- xc(view(xyz, :, 3))
    # dy = view(xyz, :, 2) .- yc(view(xyz, :, 3))
	dxy = hcat(([i-k*c1+x0, j-k*c2+y0] for (i,j,k) in eachcol(xyz))...)
	dx = view(dxy, 1, :)
    dy = view(dxy, 2, :)
    return bg .- A * exp.( -(1 .-sqrt.(( dx .* ct .+ dy .* st ).^2/rx^2+(dx .* st .- dy .* ct).^2/ry^2)).^2/σ^2)
end

function fitellipse3d(imgs::AbstractArray, mask::BitArray, edge::Matrix{Int};
                     verbose=true, keepinitialonerror=true, secondfit=true)
    r,c,h = size(imgs)
    # find ecllipse in edge
    coords = hcat(([i.I...] for i in findall(edge.>0))...)
    E1 = fitellipsedirect(coords) |> canonical
    verbose && @debug "Initial fit" E1

    # get all points within the outer boundry
    idxs = findall(mask)
    coords2 = hcat(([ci.I...] for ci in idxs)...)
    z = @view imgs[idxs]

    # refine ellipse paramaters using outer segment points
    # and edge ellipse estimate
    mmz = extrema(z)
    p0 = prepare_initial_point([E1[1:5]..., π/40, mmz..., 0.01])
    fit = curve_fit(gaussellipse3d, coords2, z, p0, autodiff=:forwarddiff)
    E = fit.param
    verbose && @debug "LSQ fit 1" p0 E

    # refine again
    lb = [c/2-3.0, r/2-3.0, 0.0, 0.0, -π, -π/2,   0,   0, 0.001]
    ub = [c/2+3.0, r/2+3.0, Inf, Inf,  π,  π/2, Inf, Inf, 1.000]
    p0 = prepare_initial_point(E)
    E = try
        if secondfit
            fit = curve_fit(gaussellipse3d, coords2, z, p0, lower=lb, upper=ub, autodiff=:forwarddiff)
            fit.param
        else
            p0
        end
    catch ex
        p0 = prepare_initial_point([E1[1:5]..., 0, mmz[1], mmz[2]/2, 0.01])
        if keepinitialonerror
            p0
        else
            @debug "Error. Trying with constraints." p0
            fit = curve_fit(gaussellipse3d, coords2, z, p0, lower=lb, upper=ub)
            fit.param
        end
    end
    verbose && secondfit && @debug "LSQ fit 2" p0 E

    return E
end

function getellipse3d(imgs::AbstractArray; verbose=true, secondfit=true)
    segs = segment3.(eachslice(imgs, dims=3))
    mask = cat(labels_map.(segs)..., dims=3) .!= 1
    fitellipse3d(imgs, mask, sum(edge3.(segs)); verbose, secondfit)
end
