using PyPlot

let
    
    k1 = 10
    k2 = 10.4

    kavg = 0.5 * (k1 + k2)
    kbeat = 0.5 * (k2 - k1) 

    xlims = pi /kbeat 

    xvals = range(-xlims, xlims, 2000)

    y = (1 .+ cos.(k1 * xvals) + cos.(k2 * xvals))
    boundingup =1 .+ 2 * cos.(kbeat .* xvals)
    boundingdown = 1 .- 2 * cos.(kbeat .* xvals)

    pygui(true)
    plot(xvals, y)
    ax = gca()

    ax[:set_xlim]([-xlims,xlims])
    plot(xvals, boundingup, color = "red", linestyle = "dashed")
    plot(xvals, boundingdown, color = "red", linestyle = "dashed")
    xlabel("Path length difference")
    ylabel("Intensity")


end