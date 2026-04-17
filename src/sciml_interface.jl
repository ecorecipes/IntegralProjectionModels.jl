"""
SciML ecosystem interface.

Converts discrete IPM problems to `DiscreteProblem` and continuous generator
problems to `ODEProblem` / `DDEProblem`.
"""

"""
    to_discrete_problem(prob::IPMProblem)

Convert an IPMProblem to a `SciMLBase.DiscreteProblem`.

The resulting problem can be solved with `OrdinaryDiffEq.FunctionMap()` or
similar discrete solvers from the SciML ecosystem.
"""
function to_discrete_problem(prob::IPMProblem{S, DensityIndependent, Deterministic}) where {S}
    K = materialize(prob.kernel)
    normalize = prob.normalize
    function ipm_step!(du, u, p, t)
        mul!(du, K, u)
        if normalize
            s = sum(du)
            if s > 0
                du ./= s
            end
        end
        return nothing
    end
    p0 = prob.p === nothing ? SciMLBase.NullParameters() : prob.p
    return SciMLBase.DiscreteProblem(ipm_step!, prob.n0, Float64.(prob.tspan), p0)
end

function to_discrete_problem(prob::IPMProblem{S, DensityDependent, Deterministic}) where {S}
    kernel_fn = prob.kernel
    normalize = prob.normalize
    function ipm_step!(du, u, p, t)
        kernel = kernel_fn(u, Int(t), p)
        K = materialize(kernel)
        mul!(du, K, u)
        if normalize
            s = sum(du)
            if s > 0
                du ./= s
            end
        end
        return nothing
    end
    p0 = prob.p === nothing ? SciMLBase.NullParameters() : prob.p
    return SciMLBase.DiscreteProblem(ipm_step!, prob.n0, Float64.(prob.tspan), p0)
end

# Wrapper to make SciML's DiscreteProblem work with in-place functions
struct DiscreteFunctionClosure{F}
    f::F
end

function (d::DiscreteFunctionClosure)(du, u, p, t)
    d.f(du, u, p, t)
end
