"""
Plot recipes for IPM visualization using RecipesBase.
"""

# Kernel heatmap
@recipe function f(::Type{Val{:kernel_heatmap}}, K::AbstractMatrix, domain::ContinuousDomain)
    z = meshpoints(domain)
    seriestype := :heatmap
    xlabel --> "Size (t)"
    ylabel --> "Size (t+1)"
    title --> "IPM Kernel"
    z, z, K
end

"""
    plot_kernel(K, domain)

Plot a materialized kernel matrix as a heatmap.
"""
@recipe function f(::Val{:ipm_kernel}, K::AbstractMatrix, domain::ContinuousDomain)
    z = meshpoints(domain)
    seriestype := :heatmap
    xlabel --> "Size (t)"
    ylabel --> "Size (t+1)"
    colorbar_title --> "Transition rate"
    z, z, K
end

# Population trajectory
@recipe function f(sol::IPMSolution)
    total_pop = [sum(u) for u in sol.u]
    xlabel --> "Time"
    ylabel --> "Total population"
    title --> "IPM Population Trajectory"
    label --> "N(t)"
    linewidth --> 2
    sol.t, total_pop
end

# Stable distribution
@recipe function f(::Val{:stable_dist}, sol::IPMSolution, domain::ContinuousDomain)
    z = meshpoints(domain)
    w = stable_distribution(sol)
    xlabel --> "Size"
    ylabel --> "Density"
    title --> "Stable size distribution"
    label --> "w(z)"
    linewidth --> 2
    fillalpha --> 0.3
    fillrange --> 0
    z, w
end
