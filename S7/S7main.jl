using CSV
using DataFrames
using PyPlot
using Peaks
using LinearAlgebra
using Statistics

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
    [24420.0, 23480.0]  #18
]

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

    data = CSV.read(file, DataFrame)
    times = data[2:end,2]
    V = data[2:end,3]

    m = zeros(length(times), 2)

    # for i = 1:length(times)
    #     m[i,1] = parse(Float64, times[i])
    #     m[i,2] = parse(Float64, V[i])
    # end
    
    m = hcat(times,V)

    return m
    
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

function stackedPlot(dataFileNumbers::Vector{Vector{Int}}, dataWnRanges::Vector{Vector{Vector{Float64}}}, correction::Function, colors::Vector{String} = color_base)

    numSubPlot = size(dataFileNumbers)[1]
    fig = figure("pyplot_subplot_touching")
    subplots_adjust(hspace=0.1)
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

        #get ridd of labels if not the last plot
        if (i != numSubPlot)
            setp(ax1.get_xticklabels(),visible=false)
        end

        #add to plots 
        numData = size(dataFileNumbers[i])[1]
        for j = 1:numData
            fileName = "S7/run" * string(dataFileNumbers[i][j]) * ".csv"
            plotData(fileName, dataWnRanges[i][j][1], dataWnRanges[i][j][2], correction, colors[colorIndex])
            if (colorIndex == length(colors))
                colorIndex = 1
            else
                colorIndex = colorIndex + 1
            end
        end

        oldaxis = ax

    end


end

function plotData(dataFileName::String, wnstart::Float64, wnstop::Float64, correction::Function, c::String)

    #dataFileName: File path name to data file
    #wnstart: monochrometer reading associated with first data points
    #wnstop: monochromator reading associated with second data points
    #correction: function that turns estimated wavenumber from monochrometer readings to true wavenumbers

    ######parameters#######
    dataSmoothingNumber = 200
    #######################

    data = getData(dataFileName)
    smoothed = smooth(data, dataSmoothingNumber)
    #normalizeData!(smoothed)

    rescaledData = rescaleX_endpoints(smoothed, wnstart, wnstop)
    correctedScaleData = hcat(correction.(rescaledData[:,1]), rescaledData[:,2])

    #plotting!
    wnmax = maximum(correctedScaleData[:,1])
    wnmin = minimum(correctedScaleData[:,1])
    pygui(true)
    ax = gca()
    ax[:set_xlim]([wnmin,wnmax])
    ax[:set_ylim]([0, maximum(smoothed[:,2])])
    plot(correctedScaleData[:,1], correctedScaleData[:,2], color = c)
    xlabel("Wavenumber " * L"($cm^{-1}$)")
    ylabel("Intensity (a.u.)")
    
end

function wavenumberCorrection(wnmono::Float64)

    #wnmono: monochromator readin of wavemnumber (cm-1)
    #retunrs: true wavenumber determined by calibration

    #current conversion function is linear
    ###### constants from linear fit on select peaks from neon spectra in run9.csv
    m =  1.0016
    b = 81.188
    ######

    return wnmono * m + b

end
