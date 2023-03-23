include("S7main.jl")


let 

    # plotData("S7/run13.csv", 23450.0, 23270.0, wavenumberCorrection, "blue")
    # plotData("S7/run11.csv", 23475.0, 23260.0, wavenumberCorrection, "#b31b1b")
    # plotData("S7/run16.csv", 24480.0, 22990.0, wavenumberCorrection, "blue")
    # plotData("S7/run15.csv", 24450.0, 23050.0, wavenumberCorrection, "#b31b1b")
    # plotData("S7/run17.csv", 24410.0, 23500.0, wavenumberCorrection, "green")
    # plotData("S7/run18.csv", 24420.0, 23480.0, wavenumberCorrection, "#b31b1b")
    # legend(["Perpendicular Polarization", "Parallel Polarization"])

    # plotData("S7/run10.csv", 25490.0, 22300.0, wavenumberCorrection, "black")
    # plotData("S7/run8.csv", 25130.0, 22330.0, wavenumberCorrection, "#b31b1b")
    # plotData("S7/run7.csv", 24800.0, 22400.0, wavenumberCorrection, "blue")
    # plotData("S7/run6.csv", 25440.0, 22380.0, wavenumberCorrection, "green")
    #legend(["Vial 1", "Vial 2", "Vial 3", "Vial 4"])


    #stackedPlot([[10], [7], [8], [6]], wavenumberCorrection, ramanScale = 0, xLim = [-1100.0, 500.0])
    #stackedPlot([[10], [7], [8], [6], [9]], x->x, ramanScale = 0)
    stackedPlot([[17,18], [21, 22]], x->x, ramanScale = 0)

    #plotting neon ones
   # stackedPlot([[9,28]], x->x)
    
    #looking at Rayleigh  peaks alone --> conclusion is to keep at 24400 +/- 10
    #need to comment out yaxis lim in plotData before running this
    #stackedPlot([[23,25,26,27]], wavenumberCorrection)



end