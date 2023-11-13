# ==============================================================================
# This file is part of the TPSA (Truncated Power Series Algebra) Julia package.
#
# Author: Jinyu Wan
# Email: wan@frib.msu.edu
# Version: 1.0
# Created Date: 11-01-2023
# Modified Date: 11-06-2023

function factorial(n::Int)
    if n < 0
        return 0
    end
    return [1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600][n + 1]
end

function doublefactorial(n::Int)
    if n < 0
        return 0
    end
    return [1, 1, 2, 3, 8, 15, 48, 105, 384, 945, 3840, 10395, 46080, 135135, 645120, 2027025, 10321920, 34459425, 185794560][n + 1]
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
