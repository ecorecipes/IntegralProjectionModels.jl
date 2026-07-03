"""
Age×size model expansion.

Expands age-generic kernel definitions into age-specific block-structured
kernels, mirroring ipmr's internal-age_size.R functionality.
"""

"""
    AgeStructure

Defines age classes for an age×size IPM.

# Fields
- `max_age`: maximum age (individuals at max_age remain there)
"""
struct AgeStructure
    max_age::Int

    function AgeStructure(max_age::Int)
        max_age > 0 || throw(ArgumentError("max_age must be positive"))
        new(max_age)
    end
end

ages(a::AgeStructure) = 1:a.max_age
n_ages(a::AgeStructure) = a.max_age

"""
    expand_age_kernels(p_func, f_func, age_structure, domain)

Expand age-generic P and F kernel functions into an age-structured MegaKernel.

# Arguments
- `p_func`: function `(age) -> AbstractSubKernel` — survival-growth kernel for each age
- `f_func`: function `(age) -> AbstractSubKernel` — fecundity kernel for each age
- `age_structure`: AgeStructure specifying age classes
- `domain`: ContinuousDomain for the size variable

Returns a block matrix where:
- P kernels appear on the sub-diagonal (aging: age a → age a+1)
- P_max_age on the diagonal for the last age class (stay at max age)
- F kernels in the first row block (all ages produce age-1 offspring)
"""
function expand_age_kernels(p_func, f_func, age_structure::AgeStructure,
        domain::ContinuousDomain)
    n_age = n_ages(age_structure)
    m = n_states(domain)
    total = m * n_age
    K = zeros(Float64, total, total)

    for a in 1:n_age
        # Block indices
        row_start = (a - 1) * m + 1
        row_end = a * m

        # Fecundity: age a produces offspring in age class 1
        F_a = materialize(f_func(a))
        K[1:m, row_start:row_end] .+= F_a

        # Survival-growth
        P_a = materialize(p_func(a))
        if a < n_age
            # Age a → age a+1 (sub-diagonal block)
            next_row_start = a * m + 1
            next_row_end = (a + 1) * m
            K[next_row_start:next_row_end, row_start:row_end] .= P_a
        else
            # Max age stays at max age (diagonal block)
            K[row_start:row_end, row_start:row_end] .+= P_a
        end
    end

    return K
end
