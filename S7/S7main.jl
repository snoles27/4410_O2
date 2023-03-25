using CSV
using DataFrames
using PyPlot
using Peaks
using LinearAlgebra
using Statistics
using DataFitting
using LsqFit

##color pallet##
color_base = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#B31B1B", "000000"]

wavenumberData = [
    [0.0,0.0],          #1
    [0.0,0.0],          #2
    [0.0,0.0],          #3
    [21770.0, 22900.0], #4
    [24100.0, 20690.0], #5
    [25440.0, 22380.0], #6
    [24800.0, 22400.0], #7
    [25130.0, 22330.0], #8
    [24550.0, 22260.0], #9
    [25490.0, 22300.0], #10
    [23475.0, 23260.0], #11
    [0.0,0.0],          #12
    [23450.0, 23270.0], #13
    [0.0, 0.0],         #14
    [24450.0, 23050.0], #15
    [24480.0, 22990.0], #16
    [24410.0, 23500.0], #17
    [24420.0, 23480.0], #18
    [24500.0, 23480.0], #19
    [24390.0, 23490.0], #20
    [24400.0, 23200.0], #21
    [24460.0, 23200.0], #22
    [24450.0, 24280.0], #23
    [0.0, 0.0],         #24
    [24470.0, 24270.0], #25
    [24450.0, 24260.0], #26
    [24450.0, 24240.0], #27
    [24100.0, 22000.0], #28
    [24099.0, 22275.0]  #29


]

labels = [
    "1",     #1
    "2",     #2
    "3",     #3
    "4",     #4
    "5",     #5
    "C6H6",  #6 
    "CHCl3", #7
    "CCl4",  #8
    "Neon 1",#9
    "CS2",   #10
    "11",    #11
    "12",    #12
    "13",    #13
    "14",    #14
    "15",    #15
    "16",    #16
    "CHCl3 (Perpendicular)", #17
    "CHCl3 (Parallel)",      #18
    "CHCl3 (Parallel)",      #19
    "CHCl3 (Perpendicular)", #20
    "CCl4 (Perpendicular)",  #21
    "CCl4 (Parallel)",       #22
    "CCl4 Rayleigh", #23
    "bad", #24
    "CHCl3 Rayleigh 1", #25
    "CHCl3 Rayleigh 2", #26
    "C6H6  Rayleigh",   #27
    "Neon 2", #28
    "Neon 3" #29
]

globe_rayleighCenter = 24440.0

function pullNumber(input)

    #parses string input for a float64

    if(input[end] == '*' || input[end] == 'f')
        return parse(Float64, input[1:end-1])
    else
        return parse(Float64, input[1:end])
    end

end

function getNeonData(file::String)

    #pulls wavenumbers and relative intensity from file

    data = CSV.read(file, DataFrame)
    wns = data[:,2]
    inten = data[:,3]
    numInten = zeros(length(inten), 1)

    for i = 1:length(inten)
        numInten[i] = pullNumber(inten[i])
    end

    return hcat(wns, numInten)

end

function getData(file::String)

    #pulls test data from "file"
    #returns [time, voltage]

    # data = CSV.read(file, DataFrame)
    # times = data[2:end,2]
    # V = data[2:end,3]

    # m = zeros(length(times), 2)

    # # for i = 1:length(times)
    # #     m[i,1] = parse(Float64, times[i])
    # #     m[i,2] = parse(Float64, V[i])
    # # end
    
    # m = hcat(times,V)

    return getData(file, 2, 3)
    
end

function getData(file::String, xindex::Int, yindex::Int)

    data = CSV.read(file, DataFrame)
    xdata = data[2:end, xindex]
    ydata = data[2:end, yindex]

    return hcat(xdata, ydata)

end

function relIntenCutOff(relInten::Float64, data::Matrix{Float64})

    #relInten: lower limit to values to be included
    #data: nx2 matrix containing [time, intensity/voltage]
    #returns: mx2 matrix where m is the number of data points above the relInensity cutoff

    m = Array{Float64}(undef, 0, 2)
    compare = maximum(data[:,2]) * relInten

    for i = 1:size(data)[1]
        if data[i,2] >= compare
            m = vcat(m, data[i,:]')
        end
    end
    return m
    
end

function rangeCutOff(range::Vector{Float64}, data::Matrix{Float64})

    #range: 2 element vecotr containing [lowerlim, upperlim]
    #data: nx2 matrix with time data in first column and value data in second column.
    #returns: rows of "data" with times within the range

    m = Array{Float64}(undef, 0, 2)

    for value in eachrow(data)
        if value[1] >= range[1] && value[1] <= range[2]
            m = vcat(m, value')
        end
    end

    return m
end

function wlen2wnum(wavelength::Float64)
    #wavelength: Wavelength in nm
    #returns: wavenumber in cm^-1

    return 1/wavelength * 10^7

end

function wnum2wlen(wavenumber::Float64)
    #wavenumber: Wavenumber in cm^-1
    #returns: Wavelength in nm

    return 1/wavenumber * 10^7

end

function smooth(data::Matrix{Float64}, averageNumber::Int)

    #performs moving average of length "averageNumber" over data in "data"

    dataNum = size(data)[1]
    smoothed = zeros(dataNum,1)
    lowEnd = round(Int, averageNumber/2)
    highEnd = averageNumber - lowEnd

    for i = 1:dataNum

        if i < lowEnd + 1
            bottom = 1
            top = i + highEnd
        elseif i > dataNum - highEnd
            top = dataNum
            bottom = i - lowEnd
        else
            bottom = i - lowEnd
            top = i + highEnd
        end

        smoothed[i] = mean(data[bottom:top, 2])
    end
        
    return hcat(data[:,1], smoothed)

end

function getPeaks(data::Matrix{Float64})

    #finds all maxima within "data"

    peakIndexs = argmaxima(data[:,2])

    peaks = zeros(length(peakIndexs),2)

    for (i, index) in enumerate(peakIndexs)
        peaks[i,:] = data[index, :]
    end
    
    return peaks

end

function rescaleX(data::Matrix{Float64}, m::Float64, b::Float64)

    #linear rescaling of time data in first column of "data" following 
    #newValue = m * oldValue + b

    rescaled = b .+ m * data[:,1]

    return hcat(rescaled, data[:,2])

end

function rescaleX_endpoints(data::Matrix{Float64,}, start::Float64, stop::Float64)

    #time recording start and end
    tstart = data[1,1]
    tend = data[end,1]

    #get line parameters for monochromator readings
    mb = linFit([tstart, tend], [start, stop], ones(2))

    return rescaleX(data, mb[1], mb[2])
end


function linFit(times::Vector{Float64}, wavenumbers::Vector{Float64}, weights::Vector{Float64})

    #performs weighted linear fit for the times wavenumber relationship
    #all input vectors must be nx1
    #returns: x = [m,b]

    W = diagm(weights)
    A = hcat(times, ones(length(times)))
    
    x = (A' * W * A)\(A' * W * wavenumbers)

    return x

end

function normalizeData!(data::Matrix{Float64})

    #rescales second column of "data" so max value is 1

    data[:,2] = data[:,2] ./ maximum(data[:,2])

end

function plotData(dataFileName::String, wnstart::Float64, wnstop::Float64, color::String)
    plotData(dataFileName, wnstart, wnstop, x -> x, color)
end

function stackedPlot(dataFileNumbers::Vector{Vector{Int}}, correction::Function; colors::Vector{String} = color_base, yaxisScaling::Float64 = 1.05, ramanScale::Int = 0, xLim::Vector{Float64} = [0.0,0.0])

    pygui(true)
    numSubPlot = size(dataFileNumbers)[1]
    fig = figure("stackedPlots", figsize=(10,10))
    subplots_adjust(hspace=0.15)
    oldaxis = gca()
    colorIndex = 1

    for i = 1:numSubPlot

        #set up subplot to share axis
        if(i == 1)
            subplot(numSubPlot * 100 + 10 + i)
        else
            subplot(numSubPlot * 100 + 10 + i, sharex=oldaxis)
        end

        #get current axis
        ax = gca()

        #get rid of labels if not the last plot
        if (i != numSubPlot)
            setp(ax.get_xticklabels(),visible=false)
        end

        #add to plots 
        numData = size(dataFileNumbers[i])[1]
        for j = 1:numData
            fileName = "S7/run" * string(dataFileNumbers[i][j]) * ".csv"
            plotData(fileName, wavenumberData[dataFileNumbers[i][j]][1], wavenumberData[dataFileNumbers[i][j]][2], correction, colors[colorIndex], yaxisScaling, ramanScale = ramanScale)
            if (colorIndex == length(colors))
                colorIndex = 1
            else
                colorIndex = colorIndex + 1
            end
        end

        if (i != numSubPlot)
            xlabel("")
        end

        #display(labels[dataFileNumbers[i]])
        legend(labels[dataFileNumbers[i]])

        oldaxis = ax

    end

    #old code for maybe setting all x scales after the fact
    # wnmax = 0.0
    # wnmin = Inf
    # for vec in dataFileNumbers
    #     for n in vec
    #         if maximum(wavenumberData[n]) > wnmax
    #             wnmax = maximum(wavenumberData[n])
    #         end
    #         if minimum(wavenumberData[n]) < wnmin
    #             wnmin = minimum(wavenumberData[n])
    #         end
    #     end
    # end
    #oldaxis[:set_xlim]([wnmin,wnmax])

    if(xLim[1] != 0.0)
        oldaxis[:set_xlim](xLim)
    end
    

    fig.canvas.draw()

end

function maxBelowRayleigh(data::Matrix{Float64}, thresh::Float64 = 0.5)

    peaks = getPeaks(data)

    justMag = peaks[:,2]
    sorted = sort(justMag, rev = true)

    return sorted[findfirst(x -> x < sorted[1] * thresh, sorted)]

end

function plotData(dataFileName::String, wnstart::Float64, wnstop::Float64, correction::Function, c::String, yaxisScaling::Float64 = 1.05; independentXScale::Bool=true, ramanScale::Int = 0)

    #dataFileName: File path name to data file
    #wnstart: monochrometer reading associated with first data points
    #wnstop: monochromator reading associated with second data points
    #correction: function that turns estimated wavenumber from monochrometer readings to true wavenumbers
    #c: Color for plot
    #yaxisScaling: number to scale y axis by 
    #independentXScale
    #raman scale: 0 --> just wavenumber, 1 --> find rayleigh peak canter and set that to 0, 2--> set center based on hardcoded value "globe_rayleighCenter"

    ######parameters#######
    dataSmoothingNumber = 50
    #######################

    # data = getData(dataFileName)
    # smoothed = smooth(data, dataSmoothingNumber)

    # rescaledData = rescaleX_endpoints(smoothed, wnstart, wnstop)
    # correctedScaleData = hcat(correction.(rescaledData[:,1]), rescaledData[:,2])

    correctedScaleData = getCorrectedData(dataFileName, wnstart, wnstop, correction, ramanScale, dataSmoothingNumber)

    #plotylim = maxBelowRayleigh(correctedScaleData)
    
    plotylim = 15.0 #override yaxis limit if neededd for specific plot

    #plotting!
    PyPlot.matplotlib[:rc]("text", usetex=false) # allow tex rendering
    rc("font", family="times", weight="normal", size = "16")
    wnmax = maximum(correctedScaleData[:,1])
    wnmin = minimum(correctedScaleData[:,1])
    pygui(true)
    ax = gca()
    if(independentXScale)
        ax[:set_xlim]([wnmin,wnmax])
    end
    ax[:set_ylim]([0, plotylim * yaxisScaling])
    plot(correctedScaleData[:,1], correctedScaleData[:,2], color = c)
    if(ramanScale != 0)
        #xlabel("Raman Shift " * "(\$ \\mathrm{cm}^{-1}\$)")
        xlabel("Raman Shift " * "(cm \$^{-1}\$)")
    else
        xlabel("Wavenumber " * "(cm \$^{-1}\$)")
    end
    ylabel("Intensity (a.u.)")
    
end

function wavenumberCorrection(wnmono::Float64)

    #wnmono: monochromator readin of wavemnumber (cm-1)
    #retunrs: true wavenumber determined by calibration

    #current conversion function is linear
    ###### constants from linear fit on select peaks from neon spectra in run9.csv
    m =  1.00105
    b = 94.94
    ######

    return wnmono * m + b

end

function findRayleighCenter(data::Matrix{Float64}, thresh::Float64 = 0.9)

    peaks = getPeaks(data)
    onlyRayleigh = Array{Float64}(undef, 0, 2)
    max = maximum(peaks[:,2])

    for row in eachrow(peaks)
        if row[2] > thresh * max
            onlyRayleigh = vcat(onlyRayleigh, row')
        end
    end

    return (maximum(onlyRayleigh[:,1]) + minimum(onlyRayleigh[:,1]))/2

end

function centerOnRayleigh(data::Matrix)

    center = findRayleighCenter(data)
    print("Rayleigh Center: ")
    println(center)
    new = centerOnRayleigh(data, center)
    # new = zeros(length(data[:,1]), length(data[1,:]))
    # new[:,1] = data[:,1] .- center
    # new[:,2] = data[:,2]

    return new

end

function centerOnRayleigh(data::Matrix, center::Float64)

    new = zeros(length(data[:,1]), length(data[1,:]))
    new[:,1] = data[:,1] .- center
    new[:,2] = data[:,2]
    
    return new

end

function getCorrectedData(fileName::String, wnstart::Float64, wnstop::Float64, correction::Function, ramanScale::Int = 0, dataSmoothingNumber::Int = 50)

    data = getData(fileName)
    smoothed = smooth(data, dataSmoothingNumber)

    rescaledData = rescaleX_endpoints(smoothed, wnstart, wnstop)
    correctedScaleData = hcat(correction.(rescaledData[:,1]), rescaledData[:,2])

    if(ramanScale == 1)
        correctedScaleData = centerOnRayleigh(correctedScaleData)
    end

    if(ramanScale == 2)
        correctedScaleData = centerOnRayleigh(correctedScaleData, globe_rayleighCenter)
    end

    return correctedScaleData

end

function getCorrectedData(dataNumber::Int, correction::Function; ramanScale::Int = 0, dataSmoothingNumber::Int = 50)

    fileName = "S7/run" * string(dataNumber) * ".csv"
    wnstart = wavenumberData[dataNumber][1]
    wnstop = wavenumberData[dataNumber][2]

    data = getData(fileName)
    smoothed = smooth(data, dataSmoothingNumber)

    rescaledData = rescaleX_endpoints(smoothed, wnstart, wnstop)
    correctedScaleData = hcat(correction.(rescaledData[:,1]), rescaledData[:,2])

    if(ramanScale == 1)
        correctedScaleData = centerOnRayleigh(correctedScaleData)
    end

    if(ramanScale == 2)
        correctedScaleData = centerOnRayleigh(correctedScaleData, globe_rayleighCenter)
    end

    return correctedScaleData

end

function trapInt(data::Matrix{Float64})

    sum = 0
    length = size(data[:,1])[1]
    for i = 1:length-1
        sum = sum + 0.5 * (data[i,2] + data[i+1,2]) * (data[i+1,1] - data[i,1])
    end
    return sum 
    
end