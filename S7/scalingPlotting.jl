include("S7main.jl")

let 

    ###########parameters############

        neon_fileName = "S7/Ne_1_300_700.csv"
        data_fileName = "S7/run9.csv"

        dataSmoothingNumber = 30

        #numbers from monochromator correspondding to the file above
        wnstart = 24550.0 #cm-1
        wnstop = 22260.0 #cm-1

        relDataPeakIntensity = 0.1

        #run6
        # peak_timerange_low = 3500.0
        # peak_timerange_high = 3600.0

        #run7
        # peak_timerange_low = 2700.0
        # peak_timerange_high = 2800.0

        #run8
        # peak_timerange_low = 3100.0
        # peak_timerange_high = 3200.0

        reference_peakWns = [22605.2, 22593.5] #data pulled from "S7/Ne_1_300_700.csv"

        reference_rayleigh = 24691

        spectralLineScale = 20

    ####endparameters##########

    #retreive neon data and remove low intensity wn 
    neondata = getNeonData(neon_fileName)
    shortenedNeon = relIntenCutOff(0.007, neondata)
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
    #cuttoff peaks by range
    # peaksCuttoff = rangeCutOff([peak_timerange_low, peak_timerange_high], peaks)
    # print("Peaks Included in fit: ") 
    # display(peaksCuttoff)

    # peakTimes = peaksCuttoff[:,1]
    # peakValues = peaksCuttoff[:,2]
    # fitLine = linFit(peakTimes, reference_peakWns, peakValues)

    # scaled_Data = rescaleX(smoothed, fitLine[1], fitLine[2])
    # scaled_peaks = rescaleX(peaks, fitLine[1], fitLine[2])
    scaled_Data = rescaleX(smoothed, roughLine[1], roughLine[2])
    scaled_peaks = rescaleX(peaks, roughLine[1], roughLine[2])

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
    xlabel("Rough Wavenumber (cm^-1)")
    ylabel("Relative Intensity")

    for row in eachrow(shortenedNeon)
        vlines(row[1], 0, row[2] * spectralLineScale, color = "red")
    end

    # plot(data4[:,1], data4[:,2])
    # plot(smoothed[:,1], smoothed[:,2])
    # scatter(peaksCuttoff[:,1], peaksCuttoff[:,2], c = "red")


end