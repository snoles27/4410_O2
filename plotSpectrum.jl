using Plots

function plotSpectrum(λmean::Float64, Δλ::Float64, N::Int)

    #Overlays a set of sinusoides with variables wavelengths

    numPoints = 10000
    xlim = λmean * (λmean + Δλ) / Δλ #beat wavelength calculated with closest wavelength to mean wavelength
    #xlim = λmean * (λmean + round(N/2) * Δλ) / (round(N/2) * Δλ) #beat wavelength calculated with the middle of the wavelength range
    #xlim = 40 * λmean
    xrange = collect(range(-xlim, xlim, length = numPoints))
    result = zeros(numPoints)

    for n = -N:N
        a = pi/(λmean + n * Δλ) #wavenumber of wavelength being used
        result += cos.(a .* xrange)
    end

    result = result./(2*N) #normalize by number of wavelengths included so average intensity is always .5
    plot(xrange, result)

end

function plotSpectrumBandWidth(λmean::Float64, bandwidth::Float64, N::Int)


    Δλ = bandwidth/N
    numPoints = 50000
    #xlim = λmean * (λmean + Δλ) / Δλ #beat wavelength calculated with closest wavelength to mean wavelength
    #xlim = λmean * (λmean + round(N/2) * Δλ) / (round(N/2) * Δλ) #beat wavelength calculated with the middle of the wavelength range
    xlim = 40 * λmean
    xrange = collect(range(-xlim, xlim, length = numPoints))
    result = zeros(numPoints)

    for n = -N:N
        a = pi/(λmean + n * Δλ) #wavenumber of wavelength being used
        result += sin.(a .* xrange)
    end

    result = result./(2*N)  #normalize by number of wavelengths included so average intensity is always .5
    plot(xrange, result, xlabel = "Light Path Difference (nm)", ylabel = "Intensity")


end


let 

    plotSpectrumBandWidth(542.0, 320.0, 1000)
    
end