using Plots
using Enzyme

x = [0.0, 0.6, 1.2, 1.3, 1.8, 1.9, 2.4]
y = sin.(x)
yd = cos.(x)

function Lan_interpolation(x0, x, y)
    n = length(x)
    if x0 < x[1] || x0 > x[n]
        error("x0 is out of range")
    end
    Px = 0.0
    for i in 1:n
        pj = 1.0
        for j in 1:n
            if j != i
                pj *= (x0 - x[j]) / (x[i] - x[j])
            end
        end
        Px += y[i] * pj
    end
    return Px
end

x_new = [0.55, 1.12, 1.27, 2.03, 2.15]
y_new = [Lan_interpolation(x_new[i], x, y) for i in 1:length(x_new)]

plot(x, y, label="original")
scatter!(x_new, y_new, label="interpolated")



