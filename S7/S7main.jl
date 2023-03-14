using CSV
using DataFrames
using PyPlot
using Peaks
using LinearAlgebra
using Statistics


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
    ylabel("Relative Intensity")
    
end

function wavenumberCorrection(wnmono::Float64)

    ###### constants from linear fit on select peaks on run9.csv
    m =  1.0016
    b = 81.188
    ######

    return wnmono * m + b

end
