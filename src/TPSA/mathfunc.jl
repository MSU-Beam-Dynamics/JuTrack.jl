# ==============================================================================
# This file is part of the TPSA (Truncated Power Series Algebra) Julia package.
#
# Author: Jinyu Wan
# Email: wan@frib.msu.edu
# Version: 1.0
# Created Date: 11-01-2023
# Modified Date: 07-02-2024

function factorial_double(n::Int)
    # factorial function for double precision
    if n <= 1
        return 1.0
    else
        return n * factorial_double(n - 1)
    end
end

function doublefactorial_double(n::Int)
    # double factorial function for double precision
    if n <= 1
        return 1.0
    else
        return n * doublefactorial_double(n - 2)
    end
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