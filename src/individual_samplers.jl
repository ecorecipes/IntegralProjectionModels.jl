"""
Analytic individual-level samplers for the parametric vital-rate types, extending
the `StructuredPopulationCore` sampler interface used by individual-based
(agent / ECS) realizations. (`sample_survives` uses the Core default, since
survival types are probability callables.)
"""

import StructuredPopulationCore: sample_growth, expected_offspring, sample_recruit

# --- Growth (sample the next trait of a survivor) ---

sample_growth(rng::AbstractRNG, g::NormalGrowth, z, domain::ContinuousDomain) =
    g.intercept + g.slope * z + g.sigma * randn(rng)

sample_growth(rng::AbstractRNG, g::LogNormalGrowth, z, domain::ContinuousDomain) =
    exp(g.intercept + g.slope * z + g.sigma * randn(rng))

# --- Expected per-parent offspring count ---

expected_offspring(f::FecundityRate, z, domain::ContinuousDomain) =
    f.establishment_prob * exp(f.intercept + f.slope * z)

expected_offspring(f::LogisticFecundityRate, z, domain::ContinuousDomain) =
    logistic(f.repr_int + f.repr_slope * z) * exp(f.seed_int + f.seed_slope * z)

# --- Recruit trait (offspring size) ---

sample_recruit(rng::AbstractRNG, f::FecundityRate, z, domain::ContinuousDomain) =
    f.recruit_mean + f.recruit_sd * randn(rng)

sample_recruit(rng::AbstractRNG, f::LogisticFecundityRate, z, domain::ContinuousDomain) =
    f.recruit_mean + f.recruit_sd * randn(rng)

sample_recruit(rng::AbstractRNG, r::RecruitmentDistribution, z, domain::ContinuousDomain) =
    r.mean + r.sd * randn(rng)
