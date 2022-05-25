using Pkg
Pkg.activate(".")

include("./utils.jl")

using Gtk.ShortNames
using ImageFiltering
using ImageView
using TestImages
using Plots
using .Utils

function wait_for_window(guidict)
    #If we are not in a REPL
    if (!isinteractive())

        # Create a condition object
        c = Condition()

        # Get the window
        win = guidict["gui"]["window"]

        # Notify the condition object when the window closes
        signal_connect(win, :destroy) do widget
            notify(c)
        end

        # Wait for the notification before proceeding ...
        wait(c)
    end
end

function lgn(stimulus;  sig_ctr=0.5, sig_srd=1.5, k_srd=0.9, alpha_mask=0.5, Vmax=275, V0=-4,
                        sig_sf=1.4, sig_u=0.3, sig_d=0.5, k_d=0.5, c50=0.1
    )
    
    # filter stimulus with receptive field
    L = customDoG(sig_ctr, sig_srd, k_srd)
    center = max.(0, L)
    sourround = max.(0, -L)

    # filter stimulus with supressive field
    H = customDoG(sig_u, sig_d, k_d)
    stimulus_sup = imfilter(stimulus, reflect(H)).^2
    
    Gsf = Kernel.gaussian((sig_sf, sig_sf), size(stimulus))
    Gsf = parent(Gsf)
    
    clocal = sqrt(sum(stimulus_sup.*Gsf))
    
    V = Vmax .* L./(c50+clocal)

    return max.(0, V.-V0)


end

stim = Utils.stimulus()
print("size of stimulus: ", size(stim), "\n")
resp = lgn(stim)
gobj = imshow(resp)

wait_for_window(gobj)
gobj = imshow(stim)

