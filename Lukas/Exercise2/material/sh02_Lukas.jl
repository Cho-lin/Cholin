Pkg.activate(".")
include("utils.jl")

using ImageFiltering
using ImageView
using Plots
using TestImages
using LinearAlgebra


function lgn(
    stimulus;
    σ_c = 0.5,
    σ_s = 1.5,
    k = 0.9,
    σ_sf = 1.4,
    σ_u = 0.3,
    σ_d = 0.5,
    k_d = 0.5,
    c50 = 0.1,
    Vmax = 273.0,
    V_0 = -6,
)
    filter_cs = Main.Utils.customDoG(σ_c, σ_s, k)
    L = imfilter(stimulus, reflect(filter_cs))
    # L_norm = Main.Utils.normalize(L)

    H = Main.Utils.customDoG(σ_u, σ_d, k_d)

    #Gsf = Kernel.gaussian((σ_sf, σ_sf), size(stimulus))
    Gsf = Kernel.gaussian(σ_sf)

    # Calculate c_local
    f = buf->sqrt(sum(dot((imfilter(buf, reflect(H)) .^ 2), parent(Gsf))))
    c_local = Main.Utils.normalize(mapwindow(f, stimulus, size(Gsf)))

    V = Vmax * (L ./ (c50 .+ c_local))

    return max.(V .- V_0, 0)
end


stimulus = Main.Utils.stimulus()
# stimulus = testimage("cameraman")
p_in = plot(Gray.(stimulus), title = "Stimulus")

out = lgn(stimulus)
p_out = plot(Gray.(Main.Utils.normalize(out)), title = "Response")

plot(p_in, p_out, layout = (1, 2))


# Plot response strength (varying radius)
rads = [0.5, 1, 2, 4, 8, 16]
outputs = Array{Float64}(undef, length(rads))

for i in 1:length(rads)
    stimulus = Main.Utils.stimulus(radius=rads[i]*Main.Utils.ppd)
    out = lgn(stimulus, Vmax=128.0, V_0=-2)
    outputs[i] = maximum(out)
end

p1 = plot(collect(zip(rads, outputs)),
 title="Response Strenght (radius)", label=false)


# Plot response strength (varying contrast)
contrasts = [0.01, 0.25, 0.50, 0.75, 1.0]
contrasts = collect(range(0, 1, 26))
outputs = Array{Float64}(undef, length(contrasts))

for i in 1:length(contrasts)
    stimulus = Main.Utils.stimulus(contrast=contrasts[i])
    out = lgn(stimulus, Vmax=128.0, V_0=-2)
    outputs[i] = maximum(out)
end

p2 = plot(collect(zip(contrasts, outputs)),
 title="Response Strenght (contrast)", label=false)

plot(p1, p2, layout=(2, 1))


# Exercise 2
img = testimage("lake")
img = Gray.(img)

img_response = lgn(img)

imshow(img_response)
