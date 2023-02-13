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

    return count

end


let

    λ1 = 380
    λ2 = 700
    
    kmax = 2 * pi /λ1
    kmin = 2 * pi /λ2
    kr = kmax - kmin
    kmid = (kmax + kmin)/2

    xrange = 15000 #nm
    xvals = range(-xrange , xrange, 1000)
    A = 1
    f(x) = A * (1 + 2 * cos(kmid * x) * sin(kr/2 * x)/(kr * x))

    threashold = .95
    edge = A^2 * (1-threashold)

    values = f.(xvals).^2
    plot(xvals, values)
    hline!([A^2 - edge], color = "red")


    #solve for the number of fringes below the line
    count = 0
    mins = argminima(values)
    for i in mins

        if(values[i] < A^2 - edge)
            count = count + 1
        end

    end

    display(count)

   
end