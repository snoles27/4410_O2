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

function getBandWidths(counts, λr, minConut, maxCount)

    bandwidths = zeros(0)

    for i = 1:length(counts)

        if counts[i] <= maxCount && counts[i] >= minConut
            append!(bandwidths, λr[i])
        end

    end

    return bandwidths

end


let

    λmid = 580.0 #nm
    λlow = 50
    λhigh = 300
    λr = range(λlow, λhigh, 3000)

    setup(bandWidth) = numFringes(λmid, bandWidth, .95)

    counts = setup.(λr)

    #DATA
    # whiteCountMin = 19
    # whiteCountMax = 20

    # greenCountMin = 57
    # greenCountMax = 65

    yellowCountMin = 43
    yellowCountMax = 50


    #setting up bandwidth calcs
    #bandwidthsWhite = getBandWidths(counts, λr, whiteCountMin, whiteCountMax)
    # bandwidthsGreen = getBandWidths(counts, λr, greenCountMin, greenCountMax)
    bandwidthsYellow = getBandWidths(counts, λr, yellowCountMin, yellowCountMax)

    PyPlot.matplotlib[:rc]("text", usetex=false) # allow tex rendering
    rc("font", family="times", weight="normal", size = "12")
    pygui(true)
    plot(counts, λr, label = "\$ \\lambda_{\\mathrm{mid}} \$ = " * string(trunc(Int, λmid)) * " nm", color = "goldenrod")
    ax = gca()

    # add new limits
    ax[:set_xlim]([minimum(counts),maximum(counts)])
    ax[:set_ylim]([λlow,λhigh])
    # y1_w = [minimum(bandwidthsWhite), minimum(bandwidthsWhite)]
    # y2_w = [maximum(bandwidthsWhite), maximum(bandwidthsWhite)]
    # x_w = [minimum(counts), whiteCountMax]
    # fill_between(x_w, y2_w, y1_w, where=(y1_w .< y2_w), color = "red", alpha = 0.3)

    # y1_g = [minimum(bandwidthsGreen), minimum(bandwidthsGreen)]
    # y2_g = [maximum(bandwidthsGreen), maximum(bandwidthsGreen)]
    # x_g = [minimum(counts), greenCountMax]
    # fill_between(x_g, y2_g, y1_g, where=(y1_g .< y2_g), color = "red", alpha = 0.3)

    y1_y = [minimum(bandwidthsYellow), minimum(bandwidthsYellow)]
    y2_y = [maximum(bandwidthsYellow), maximum(bandwidthsYellow)]
    x_g = [minimum(counts), yellowCountMax]
    fill_between(x_g, y2_y, y1_y, where=(y1_y .< y2_y), color = "red", alpha = 0.3)

    ylabel("Bandwidth (nm)")
    xlabel("Fringe Count")
    # vlines([whiteCountMin,whiteCountMax], 0, λhigh, color = "red", linestyle = "dashed", label = "Fringe Count Range")
    # title("(a) White Light")
    # vlines([greenCountMin,greenCountMax], 0, λhigh, color = "red",linestyle = "dashed", label = "Fringe Count Range")
    # title("(b) Green Filter")
    vlines([yellowCountMin,yellowCountMax], 0, λhigh, color = "red",linestyle = "dashed", label = "Fringe Count Range")
    title("(c) Yellow Filter")
    
    # vlines([minimum(counts)], minimum(bandwidthsWhite), maximum(bandwidthsWhite), color = "red", linewidth = 5, label = "Bandwith Range (White Light)")
    # vlines([minimum(counts)], minimum(bandwidthsGreen), maximum(bandwidthsGreen), color = "green", linewidth = 5, label = "Bandwith Range (Green Light)")
    
    legend()

    # plotFringes(λmid, 350.0, .95)

    display(y1_y[1])
    display(y2_y[1])

    
end