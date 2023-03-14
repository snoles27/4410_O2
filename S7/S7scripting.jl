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
    # legend(["Vial 1", "Vial 2", "Vial 3", "Vial 4"])

    xdata = [1, 2, 3, 4] * 0.2
    ydata = [1, 2, 3, 4] * 0.2
    fig = figure("pyplot_subplot_touching",figsize=(10,10))
    subplots_adjust(hspace=0.1) # Set the vertical spacing between axes
    subplot(311) # Create the 1st axis of a 3x1 array of axes
    ax1 = gca()
    setp(ax1.get_xticklabels(),visible=false) # Disable x tick labels
    PyPlot.title("Title")
    ylim(0.0,1.0) # Set the y-limits from 0.0 to 1.0
    subplot(312,sharex=ax1) # Create the 2nd axis of a 3x1 array of axes
    ax2 = gca()
    setp(ax2.get_xticklabels(),visible=false) # Disable x tick labels
    plot(xdata,ydata)
    grid("on")
    ylabel("Log Scale")
    yticks(0.1:0.2:0.9)
    ylim(0.0,1.0)
    subplot(313,sharex=ax2) # Create the 3rd axis of a 3x1 array of axes
    ax3 = gca()
    plotData("S7/run18.csv", 24420.0, 23480.0, wavenumberCorrection, "#b31b1b")
    grid("on")
    xlabel("Log Scale")
    yticks(0.1:0.2:0.9)
    ylim(0.0,1.0)
    fig.canvas.draw() # Update the figure
    suptitle("3x1 Shared Axis")

end