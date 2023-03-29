
include("S7main.jl")

function linear(x::Float64, m::Float64, b::Float64)
    return m * x + b
end

function linear(x::Vector{Float64}, m::Float64, b::Float64)
    return m .* x .+ b
end

function linear(x::Vector{Float64}, mb::Vector{Float64})
    return mb[1] .* x .+ mb[2]
end

function linear(x::Float64, mb::Vector{Float64})
    return mb[1] * x + mb[2]
end

function quadratic(x::Vector{Float64}, abc::Vector{Float64})
    return abc[1] .* x.^2 + abc[2] .* x .+ abc[3]
end

function lorentzian(x::Vector{Float64}, x0gb::Vector{Float64})
    #return x0gb[3] .* x0gb[2]./((x .- x0gb[1]).^2 .+ x0gb[2]^2) .+ x0gb[4]
    #return x0gb[3] .* x0gb[2]./((x .- x0gb[1]).^2 .+ x0gb[2]^2) .+ x0gb[4]
    return x0gb[3] ./ (1 .+ ((x .- x0gb[1])./x0gb[2]).^2) .+ x0gb[4]
end

function BWF(x::Vector{Float64}, x0gIQb::Vector{Float64})

    return x0gIQb[3] .* (1 .+ (x .- x0gIQb[1])/(x0gIQb[2]*x0gIQb[4])).^2 ./ (1 .+ ((x .- x0gIQb[1])/x0gIQb[2]).^2) .+ x0gIQb[5]

end

function fitLorentzian(dataNumber::Int, range::Vector{Float64} = [0.0,0]; ramanScale::Int = 0, includePlot::Bool = false)

    longdata = getCorrectedData(dataNumber, wavenumberCorrection, ramanScale = ramanScale)
    if(range[1] == 0)
        data = longdata
    else
        data = rangeCutOff(range, longdata)
    end

    xvals = data[:,1]
    y = data[:,2]

    fit = curve_fit(lorentzian, xvals, y, [mean([xvals[1], xvals[end]]), 1.0, 1.0, 0.0])
    param = coef(fit)

    if(includePlot)
        pygui(true)
        xplot = collect(LinRange(xvals[1], xvals[end], 1000))
        plot(xvals, y)
        plot(xplot, lorentzian(xplot, param))
        xlabel("Wavenumber " * "(cm \$^{-1}\$)")
        ylabel("Intensity (a.u.)")
        vlines(param[1], 0, .23, "red")
        legend(["Data", "Lorentzian Fit", "Peak Center = " * string(round(param[1]))* "cm \$^{-1}\$"])
    end

    return param, stderror(fit) 
    
end

function fitBWF(dataNumber::Int, range::Vector{Float64}; ramanScale::Int = 0, includePlot::Bool = false)

    longdata = getCorrectedData(dataNumber, wavenumberCorrection, ramanScale = ramanScale)
    data = rangeCutOff(range, longdata)

    xvals = data[:,1]
    y = data[:,2]

    fit = curve_fit(BWF, xvals, y, [mean(range), 5.0, 1.0, 1.0, 0.0])
    param = coef(fit)

    if(includePlot)
        xplot = collect(LinRange(range[1], range[2], 1000))
        plot(xvals, y)
        plot(xplot, BWF(xplot, param))
    end

    return param, stderror(fit) 
    
end

let 

    

    ## Fitting Peaks ## 

    #C6H6 peak ranges
    #range = [23725.0, 23850.0] #peak1 run10
    #range = [23575.0, 23700.0] #peak2 run10

    ##CCl4 peak ranges
    #range = [24180.0, 24280.0] #peak1 run8
    #range = [24080.0, 24180.0] #peak2 run8
    #range = [23900.0, 24050.0] #peak3 run8
    range = [23550.0, 23750.0] #peak4 run8
    #range = [24621.0, 24694.0]  #anti 1 run8
    #range = [24720.0, 24789.0]  #anti 2 run8
    #range = [24861.0, 24956.0]  #anti 3 run 8 (CCl4)

    ##CHCl3 Peak ranges
    #range = [24140.0, 24230.0] #peak 1 run 7 (CHCl3)
    # range = [24030.0, 24120.0] # peak 2 
    # range = [24066.0, 24090.0] #right sub peak of peak 2
    # range = [24050.0,24065.0] #left sub peak of peak 2
    # range = [23729.0,23820.0]  # peak 3
    # range = [23620.0, 23730.0] #peak 4
    # range = [24670.0, 24760.0] #anti 1
    # range = [24780.0, 24860.0] #anti 2

    ##C6H6 Peak ranges
    #range = [23796.0,23866.0] #peak 1
    #range = [23391.0,23507.0] #peak 2
    #range = [23211.0,23302.0] #peak 3
    #range = [22785.0,22900.0] #peak 4

    fitLorentzian(8, range, includePlot = true)
    #fitBWF(10, range, includePlot = true)

    #fitLorentzian(27, includePlot = true)

    ###OLD###

    ### Getting averages out of run31 for repeatablilty measurements
    # rawdata = CSV.read("S7/run31.csv", DataFrame)
    # xdata = rawdata[2:end, 2]
    # ydata = rawdata[2:end, 3]
    # plot(xdata,ydata)

    # data = hcat(xdata, ydata)

    # ll1 = mean(rangeCutOff([20,61.0], data)[:,2])
    # ll2 = mean(rangeCutOff([209,293.0], data)[:,2])
    # ll3 = mean(rangeCutOff([407,497.0], data)[:,2])

    # display([ll1,ll2,ll3])

    # T1 = mean(rangeCutOff([85,167.0], data)[:,2])
    # T2 = mean(rangeCutOff([318,386.0], data)[:,2])
    # T3 = mean(rangeCutOff([519,639.0], data)[:,2])

    # display([T1,T2,T3])
    
    ### Getting averages out of run30 for repeatablilty measurements
    # rawdata = CSV.read("S7/run30.csv", DataFrame)
    # xdata = rawdata[2:end, 2]
    # ydata = rawdata[2:end, 3]

    # data = hcat(xdata, ydata)

    # ll1 = mean(rangeCutOff([0,62.0], data)[:,2])
    # ll2 = mean(rangeCutOff([312,401.0], data)[:,2])
    # ll3 = mean(rangeCutOff([551,614.0], data)[:,2])

    # display([ll1,ll2,ll3])

    # T1 = mean(rangeCutOff([142,278.0], data)[:,2])
    # T2 = mean(rangeCutOff([433,522.0], data)[:,2])
    # T3 = mean(rangeCutOff([636,740.0], data)[:,2])
   
    # display([T1,T2,T3])

    # longData = getCorrectedData(8, wavenumberCorrection, ramanScale = 0)
    # plot(longData[:,1], longData[:,2])

    #####general fitting before function writing######

    # longdata = getCorrectedData(10, wavenumberCorrection, ramanScale = 0)
    # data = rangeCutOff(range, longdata)

    # xvals = data[:,1]
    # y = data[:,2]

    # plot(xvals, y)

    # fit = curve_fit(lorentzian, xvals, y, [mean(range), 1.0, 1.0, 0.0])
    # param = coef(fit)
    # display(param)
    # display(stderror(fit))

    # xplot = collect(LinRange(range[1], range[2], 1000))

    # plot(xplot, lorentzian(xplot, param))


    # ##Getting linear fit 

    # data = getData("S7/wnrelation.csv", 1, 2)

    # xvals = data[:,1]
    # y = data[:,2]
    # pygui(true)
    # scatter(xvals, y)

    
    # fit = curve_fit(linear, xvals, y, [1.0, 40.0])
    # mb = coef(fit)
    # display(mb)
    # display(stderror(fit))

    # xplot = LinRange(minimum(xvals), maximum(xvals), 10)
    # yplot = mb[1] .* xplot .+ mb[2]

    # plot(xplot, yplot)

    # avgResid = mean(abs.(fit.resid))


    # #Old method using DataFitting --> Couldn't figure out uncertainty on the Measures object
    # paramGuess = [1.0, 0.0]
    # pygui(true)
    # scatter(xvals, y)

    # dom = Domain(xvals)
    # data = Measures(y)
    
    # testingModel = Model(:comp1 => FuncWrap(linear, paramGuess...))
    # prepare!(testingModel, dom, :comp1)
    # result1 = fit!(testingModel, data)


end