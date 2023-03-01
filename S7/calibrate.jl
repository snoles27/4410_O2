using CSV
using DataFrames
using PyPlot
using Peaks
using Statistics


function pullNumber(input)

    if(input[end] == '*')
        return parse(Float64, input[1:end-1])
    else
        return parse(Float64, input[1:end-1])
    end

end

function getNeonData(file::String)

    data = CSV.read(file, DataFrame)
    wlens = data[:,1]
    inten = data[:,4]
    numInten = zeros(length(inten), 1)

    for i = 1:length(inten)
        numInten[i] = pullNumber(inten[i])
    end

    return hcat(wlens, numInten)

end

function getData(file::String)
    data = CSV.read(file, DataFrame)
    times = data[2:end-5,2]
    V = data[2:end-5,3]

    m = zeros(length(times), 2)

    for i = 1:length(times)
        m[i,1] = parse(Float64, times[i])
        m[i,2] = parse(Float64, V[i])
    end

    return m
    
end

function relIntenCutOff(relInten::Float64, data::Matrix{Float64})

    m = Array{Float64}(undef, 0, 2)
    compare = maximum(data[:,2]) * relInten

    for i = 1:size(data)[1]
        if data[i,2] >= compare
            m = vcat(m, data[i,:]')
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

    peakIndexs = argmaxima(data[:,2])

    peaks = zeros(length(peakIndexs),2)

    for (i, index) in enumerate(peakIndexs)
        peaks[i,:] = data[index, :]
    end
    
    return peaks

end

function rescaleX(data::Matrix{Float64}, min::Float64, max::Float64)

    slope = (max - min)/(data[end,1] - data[1,1])

    rescaled = min .+ slope * data[:,1]

    return hcat(rescaled, data[:,2])

end

let 

    neondata = getNeonData("S7/Ne_I.csv")
    shortenedData = relIntenCutOff(0.1, neondata)

    data4 = getData("S7/run4.csv")
    wnmin = 21770.0 #cm-1
    wnmax = 22900.0 #cm-1

    wlenmax = wnum2wlen(wnmin)
    wlenmin = wnum2wlen(wnmax)

    smoothed = smooth(data4, 30)
    peaks = getPeaks(smoothed)
    peaksCuttoff = relIntenCutOff(0.2, peaks)

    scaleChange = rescaleX(smoothed, wlenmin, wlenmax)
    scaleChange_raw = rescaleX(data4, wlenmin, wlenmax)



    pygui(true)
    plot(scaleChange_raw[:,1], scaleChange_raw[:,2])
    plot(scaleChange[:,1], scaleChange[:,2])
    for row in eachrow(shortenedData)
        vlines(row[1], 1, 1 + row[2] * .01, color = "red")
    end





    # plot(data4[:,1], data4[:,2])
    # plot(smoothed[:,1], smoothed[:,2])
    # scatter(peaksCuttoff[:,1], peaksCuttoff[:,2], c = "red")


end