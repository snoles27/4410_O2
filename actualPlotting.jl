using PyPlot
using Peaks

function numFringes(λmid::Float64, λr::Float64, threshold::Float64)

    kmid = 2 * pi / λmid
    kr = 2 * pi * λr / (λmid^2 - .25 * λr^2)

    xrange = 20000 #nm
    xvals = range(-xrange , xrange, 2000)
    A = 1
    f(x) = A * (1 + 2 * cos(kmid * x) * sin(kr/2 * x)/(kr * x))

    edge = A^2 * (1-threshold)

    values = f.(xvals).^2
    #solve for the number of fringes below the line
    count = 0
    mins = argminima(values)
    for i in mins

        if(values[i] < A^2 - edge)
            count = count + 1
        end

    end

    # plot(xvals, values)
    # hline!([A^2 - edge], color = "red")

    return count

end

function plotFringes(λmid::Float64, λr::Float64, threshold::Float64)

    kmid = 2 * pi / λmid
    kr = 2 * pi * λr / (λmid^2 - .25 * λr^2)

    xrange = 10000 #nm
    xvals = range(-xrange , xrange, 2000)
    A = 1
    f(x) = A * (1 + 2 * cos(kmid * x) * sin(kr/2 * x)/(kr * x))

    edge = A^2 * (1-threshold)

    values = f.(xvals).^2
    #solve for the number of fringes below the line
    count = 0
    mins = argminima(values)
    storeMinLoc = zeros(0)
    storeMinValue = zeros(0)
    for i in mins

        if(values[i] < A^2 - edge)
            count = count + 1
            append!(storeMinLoc, xvals[i])
            append!(storeMinValue, values[i])
        end

    end

    pygui(true)
    plot(xvals, values)
    hlines([A^2 - edge], -xrange, xrange, colors = "red", label = "Counting Threshold")
    scatter(storeMinLoc, storeMinValue, c = "orange", s = 15, label = "Countable Fringes")
    xlabel("Pathlength Difference (nm)")
    ylabel("Relative Intensity")
    legend()



end


let

    λmid = 540.0 #nm
    λr = range(100, 600, 1000)

    setup(bandWidth) = numFringes(λmid, bandWidth, .95)

    counts = setup.(λr)

    pygui(true)
    plot(counts, λr)
    ylabel("Bandwidth (nm)")
    xlabel("Fringe Count")
    vlines([19,20], 100, 600, color = "red", label = "Observed Fringe Counts")
    legend()

    # plotFringes(λmid, 350.0, .95)

    
end