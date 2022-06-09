using Pkg
Pkg.activate(".")
include("getgabor.jl")


using ImageFiltering
using TestImages
using ImageView
using Plots
using TestImages
using LinearAlgebra
using Images


function normalize_input(input)
    return (input .- minimum(input)) ./ (maximum(input) - minimum(input))
end


#################
#               #
#  Exercise 1   #
#               #
#################
function v1(input)
    act_lgn = imfilter(input, Kernel.DoG((0.7, 2.0)))

    orientations = collect(0:15:165)
    act_v1_θ = []
    for i in 1:length(orientations)
        gabor = GetGabor.getgabor(0.1, 5, orientations[i])

        response = imfilter(Float64.(act_lgn), gabor, Algorithm.FIR())

        act_v1_θ = [act_v1_θ..., abs.(response)]
    end

    return act_v1_θ
end


img = load("/Users/lukas/Documents/Uni/Semester_2/ViMaM/Exercises/Exercise3/LiCircleStim.png")
act_v1_θ = v1(Gray.(img))


#################
#               #
# Exercise 2.1  #
#               #
#################
function longrangefilter(θ; σ=5, rmax=60, κ=50)
    s = 2*rmax + 4*σ + 1
    center = rmax + 2*σ + 1

    mask = zeros((s, s))

    for sx in collect(1:s)
        for sy in collect(1:s)
            r = hypot(sx-center, sy-center)
            φ = atan(sx-center, sy-center)

            B_ang_θ = exp(κ * cos(φ + θ))
            if abs(r) > rmax
                B_rad = exp(-((abs(r) - rmax)^2 / 2*σ^2))
            else
                B_rad = 1
            end

            mask[sx, sy] = B_ang_θ * B_rad
        end
    end

    return normalize_input(mask)
end

mask = longrangefilter(pi)
imshow(mask)
