# Author: Christian Jarvers
#
# Copyright 2020 Christian Jarvers (MIT License)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
# documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and # to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of 
# the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO 
# THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.

module Utils

using ImageFiltering

export customDoG, normalize, ppd

"""
    customDoG(σᵢₙ, σₒᵤₜ, k)

Creates a difference-of-Gaussian filter, where the central Gaussian has standard deviation `σᵢₙ`,
the surround has standard deviation `σₒᵤₜ` and is scaled by factor `k`.
"""
function customDoG(σᵢₙ, σₒᵤₜ, k)
    outer = Kernel.gaussian(σₒᵤₜ)
    inner = Kernel.gaussian((σᵢₙ, σᵢₙ), size(outer))
    return inner .- k * outer
end

"""
    normalize(img)

Normalizes `img` to be in range [0, 1].
"""
function normalize(img)
    return (img .- minimum(img)) ./ (maximum(img) - minimum(img))
end

"""
    ppd

Pixels per degree of visual input.
"""
ppd = 10.0

"""
    stimulus(; imsize, radius, contrast, θ)

Generates a visual stimulus display of size `imsize`. A central circle of radius `radius`
is filled with a sine grating with spatial frequency `θ`, with the intensity of the peaks
scaled by `contrast`.
"""
function stimulus(; imsize=Int(40*ppd), radius=20*ppd, contrast=1.0, θ=0.24/ppd)
    xs = repeat(-imsize//2:imsize//2, inner=(1, Int(imsize//2)*2+1))
    ys = xs'
    circle = sqrt.(xs.^2 + ys.^2) .≤ radius
    img = sin.(ys .* θ .* π) .* contrast .* circle
end

end # module Utils
