# ==============================================================================
# This file is part of the TPSA (Truncated Power Series Algebra) Julia package.
#
# Author: Jinyu Wan
# Email: wan@frib.msu.edu
# Version: 1.0
# Created Date: 11-01-2023
# Modified Date: 11-06-2023

# function factorial(n::Int)
#     if n < 0
#         return 0
#     end
#     return [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600][n + 1]
# end

function doublefactorial(n::Int)
    if n < 0
        return 0.0
    end
    return [1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0][n + 1]
end

function binomial(n::Int, m::Int)
    if n <= 0 || m > n || m < 0
        return 0
    end
    ml = m > n ÷ 2 ? n - m : m
    base_cases = [1, n, n * (n - 1) ÷ 2, n * (n - 1) * (n - 2) ÷ 6]
    if ml >= 0 && ml <= 3
        return base_cases[ml + 1]
    end
    if n <= 12
        return factorial(n) ÷ (factorial(m) * factorial(n - m))
    else
        return (n * binomial(n - 1, ml - 1)) ÷ ml
    end
end


function gcd(a::T, b::T) where T
    b == 0 ? a : gcd(b, a % b)
end

function fraction_reduction(a::T, b::T) where T
    if a == 0 || b == 0 
        return
    end
    temp = gcd(abs(a), abs(b))
    a /= temp
    b /= temp
end
