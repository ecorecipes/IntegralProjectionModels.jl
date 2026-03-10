"""
SciML ecosystem interface.

Converts IPMProblems to SciMLBase.DiscreteProblem for use with
DifferenceEquations.jl callbacks, ensemble simulations, etc.
"""

"""
    to_discrete_problem(prob::IPMProblem)

Convert an IPMProblem to a `SciMLBase.DiscreteProblem`.

The resulting problem can be solved with `OrdinaryDiffEq.FunctionMap()` or
similar discrete solvers from the SciML ecosystem.
"""
function to_discrete_problem(prob::IPMProblem{S, DensityIndependent, Deterministic}) where {S}
    K = materialize(prob.kernel)
    function ipm_step!(du, u, p, t)
        mul!(du, K, u)
        if prob.normalize
            s = sum(du)
            if s > 0
                du ./= s
            end
        end
    end
    return SciMLBase.DiscreteProblem(
        SciMLBase.DiscreteFunctionClosure(ipm_step!),
        prob.n0,
        Float64.(prob.tspan),
        prob.p)
end

function to_discrete_problem(prob::IPMProblem{S, DensityDependent, Deterministic}) where {S}
    function ipm_step!(du, u, p, t)
        kernel = prob.kernel(u, Int(t), p)
        K = materialize(kernel)
        mul!(du, K, u)
        if prob.normalize
            s = sum(du)
            if s > 0
                du ./= s
            end
        end
    end
    return SciMLBase.DiscreteProblem(
        SciMLBase.DiscreteFunctionClosure(ipm_step!),
        prob.n0,
        Float64.(prob.tspan),
        prob.p)
end

# Wrapper to make SciML's DiscreteProblem work with in-place functions
struct DiscreteFunctionClosure{F}
    f::F
end

function (d::DiscreteFunctionClosure)(du, u, p, t)
    d.f(du, u, p, t)
end
