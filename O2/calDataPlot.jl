using PyPlot
using CSV 
using DataFrames

let

    data = CSV.read("CalData.csv", DataFrame)

    micDistance = data[:,1]
    mirrorDistance = data[:,2]
    errorBar = .089 .* micDistance /2

    x = [.207, .0704] #line data [m,b] from google sheets fit

    xplotMin = 75
    xplotMax = 155
    linex = range(xplotMin, xplotMax, 10)

    line = x[1] * linex .+ x[2]

    PyPlot.matplotlib[:rc]("text", usetex=false) # allow tex rendering
    rc("font", family="times", weight="normal", size = "12")
    pygui(true)
    # scatter(micDistance, mirrorDistance)
    xlabel("Micrometer Distance (\$\\mu \$m)")
    ylabel("Mirror Distance (\$\\mu \$m)")
   
    errorbar(micDistance, mirrorDistance, yerr=errorBar, fmt="o",markersize=8, capsize=4)
    ax = gca()
    ax[:set_ylim]([minimum(mirrorDistance)-5 ,maximum(mirrorDistance) + 7])
    ax[:set_xlim](xplotMin, xplotMax)
    plot(linex, line, label = "Linear Fit, ratio = " * string(x[1]))
    legend()
    

end