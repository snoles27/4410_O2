using CSV
using DataFrames


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

let 

    data = getNeonData("S7/Ne_I.csv")



end