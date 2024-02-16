using Enzyme
Enzyme.API.runtimeActivity!(true)
using BenchmarkTools
function sum_mt(len)
    s = 0.0
    Threads.@threads for i in 1:10
        s += len
    end
    return s
end
function sum_no_mt(len)
    s = 0.0
    for i in 1:10
        s += len
    end
    return s
end

grad_no_mt = autodiff(Forward, sum_no_mt, Duplicated(3.0, 1.0))
println("gradient is $grad_no_mt")
grad_mt = autodiff(Forward, sum_mt, Duplicated(3.0, 1.0))
println("gradient is $grad_mt")