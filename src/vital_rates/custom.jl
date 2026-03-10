"""
    CustomVitalRate{F,P}

Wraps any callable as a vital rate. The function `func` is called as
`func(args...; params...)` or simply `func(args...)`.

# Examples
```julia
# A size-dependent vital rate
vr = CustomVitalRate((z) -> 0.5 * z, nothing)
vr(1.0)  # 0.5

# With parameters
vr = CustomVitalRate((z, p) -> p.a + p.b * z, (a=1.0, b=0.5))
vr(2.0)  # calls func(2.0, (a=1.0, b=0.5))
```
"""
struct CustomVitalRate{F, P} <: AbstractVitalRate
    func::F
    params::P
end

CustomVitalRate(func) = CustomVitalRate(func, nothing)

function (v::CustomVitalRate{F, Nothing})(args...) where {F}
    return v.func(args...)
end

function (v::CustomVitalRate{F, P})(args...) where {F, P}
    return v.func(args..., v.params)
end
