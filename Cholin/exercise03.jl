using Pkg
Pkg.activate(".")


function longrangefilter(θ, σ=5, rmax=60, κ=50)
    s = 2*rmax + 4*θ +1
    r = hypot(x,y)
    ϕ = atan(y,x)
    function Bθ(ϕ, r)

        function Bangθ(ϕ)
            return exp(κ*cos(ϕ+θ))
        end
        function Brad(r)
            if abs(r)>rmax
                return exp(-((abs(r)-rmax)^2)/(2*σ^2))
            else
                return 1
        end
        return Bangθ(ϕ)*Brad(r) #elementwise?
    end
    return #s X s filter
end

end

#generate 2 filter banks
right =
