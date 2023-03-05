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

let 

    ###########parameters############

        neon_fileName = "S7/Ne_1_300_700.csv"
        data_fileName = "S7/run7.csv"

        dataSmoothingNumber = 30

        #numbers from monochromator correspondding to the file above
        wnstart = 24800.0 #cm-1
        wnstop = 22400.0 #cm-1

        relDataPeakIntensity = 0.2

        #run6
        # peak_timerange_low = 3500.0
        # peak_timerange_high = 3600.0

        #run7
        peak_timerange_low = 2700.0
        peak_timerange_high = 2800.0

        reference_peakWns = [22605.2, 22593.5] #data pulled from "S7/Ne_1_300_700.csv"

        reference_rayleigh = 24691

        spectralLineScale = 20

    ####endparameters##########

    #retreive neon data and remove low intensity wn 
    neondata = getNeonData(neon_fileName)
    shortenedNeon = relIntenCutOff(0.01, neondata)
    normalizeData!(shortenedNeon)

    #retrieve data file
    data = getData(data_fileName)
    smoothed = smooth(data, dataSmoothingNumber)
    normalizeData!(smoothed)

    #time recording start and end
    tstart = smoothed[1,1]
    tend = smoothed[end,1]

    #get line parameters for monochromator readings
    roughLine = linFit([tstart, tend], [wnstart, wnstop], ones(2))

    #find peaks of smoothed data
    peaks = getPeaks(smoothed)
    #remove peaks below a certain relative intensity
    peaks = relIntenCutOff(relDataPeakIntensity, peaks)
    peaksCuttoff = rangeCutOff([peak_timerange_low, peak_timerange_high], peaks)
    print("Peaks Included in fit: ") 
    display(peaksCuttoff)

    peakTimes = peaksCuttoff[:,1]
    peakValues = peaksCuttoff[:,2]
    fitLine = linFit(peakTimes, reference_peakWns, peakValues)

    scaled_Data = rescaleX(smoothed, fitLine[1], fitLine[2])
    scaled_peaks = rescaleX(peaks, fitLine[1], fitLine[2])
    # scaled_Data = rescaleX(smoothed, roughLine[1], roughLine[2])
    # scaled_peaks = rescaleX(peaks, roughLine[1], roughLine[2])

    #plotting!
    wnmax = maximum(scaled_Data[:,1])
    wnmin = minimum(scaled_Data[:,1])
    pygui(true)
    ax = gca()
    ax[:set_xlim]([wnmin,wnmax])
    ax[:set_ylim]([0, 1.05])
    #plot(scaled_Data_raw[:,1], scaled_Data_raw[:,2])
    plot(scaled_Data[:,1], scaled_Data[:,2], color = "blue")
    scatter(scaled_peaks[:,1], scaled_peaks[:,2], c = "green")

    for row in eachrow(shortenedNeon)
        vlines(row[1], 0, row[2] * spectralLineScale, color = "red")
    end

    # plot(data4[:,1], data4[:,2])
    # plot(smoothed[:,1], smoothed[:,2])
    # scatter(peaksCuttoff[:,1], peaksCuttoff[:,2], c = "red")


end