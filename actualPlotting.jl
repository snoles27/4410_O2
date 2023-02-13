using Plots
using Peaks

function numFringes(λmid::Float64, λr::Float64, threshold::Float64)

    kmid = 2 * pi / λmid
    kr = 2 * pi * λr / (λmid^2 - .25 * λr^2)

    xrange = 30000 #nm
    xvals = range(-xrange , xrange, 2000)
    A = 1
    f(x) = A * (1 + 2 * cos(kmid * x) * sin(kr/2 * x)/(kr * x))

    edge = A^2 * (1-threshold)

    values = f.(xvals).^2
    #solve for the number of fringes below the line
    count = 0
    mins = argminima(values)
    for i in mins

        if(values[i] < A^2 - edge)
            count = count + 1
        end

    end

    # plot(xvals, values)
    # hline!([A^2 - edge], color = "red")

    return count

end


let

    λmid = 540.0 #nm
    λr = range(20, 500, 100)

    setup(bandWidth) = numFringes(λmid, bandWidth, .95)

    counts = setup.(λr)

    plot(λr, counts)
    
    
end