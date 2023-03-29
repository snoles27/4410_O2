include("S7main.jl")


let 
    
    ### Working on deplorization ratio ###
    #stackedPlot([[17,18], [21, 22]], wavenumberCorrection, ramanScale = 2)

    # ##CCl4 deplorization ratio calcs
    # ccl4_ll = getCorrectedData(22, wavenumberCorrection)
    # ccl4_T = getCorrectedData(21, wavenumberCorrection)

    # range_ccl4 = [[24190.0, 24250.0], [24096.0, 24155.0], [23942., 24021.], [23615., 23711.]]
    # Ill_ccl4 = zeros(4)
    # IT_ccl4 = zeros(4)
    # p = zeros(4) 

    # for i = 1:4
    #     Ill_ccl4[i] = trapInt(rangeCutOff(range_ccl4[i], ccl4_ll))
    #     IT_ccl4[i] = trapInt(rangeCutOff(range_ccl4[i], ccl4_T))
    #     p[i] = IT_ccl4[i]/Ill_ccl4[i]
    # end

    # println(p)

    ##chcl3 depolarization ratios calc
    # chcl3_ll = getCorrectedData(18, wavenumberCorrection)
    # chcl3_T = getCorrectedData(17, wavenumberCorrection)

    # range_chcl3 = [[24144.0, 24207.0], [24038.0, 24100.0], [23737., 23801.], [23636., 23713.]]
    # Ill_chcl3 = zeros(4)
    # IT_chcl3 = zeros(4)
    # p = zeros(4)

    # for i = 1:4
    #     Ill_chcl3[i] = trapInt(rangeCutOff(range_chcl3[i], chcl3_ll))
    #     IT_chcl3[i] = trapInt(rangeCutOff(range_chcl3[i], chcl3_T))
    #     p[i] = IT_chcl3[i]/Ill_chcl3[i]
    # end

    # println(p)


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


    stackedPlot([[10], [7], [8], [6]], wavenumberCorrection, ramanScale = 0)
    #stackedPlot([[10], [7], [8], [6]], wavenumberCorrection, ramanScale = 2)
    #stackedPlot([[10], [7], [8], [6], [9]], x->x, ramanScale = 0)
    #stackedPlot([[17,18], [21, 22]], wavenumberCorrection, ramanScale = 2)

    #plotting neon ones
    #stackedPlot([[9,28,29]], x->x)
    
    #looking at Rayleigh  peaks alone --> conclusion is to keep at 24400 +/- 10
    #need to comment out yaxis lim in plotData before running this
    #stackedPlot([[23,25,26,27]], wavenumberCorrection)
    
end