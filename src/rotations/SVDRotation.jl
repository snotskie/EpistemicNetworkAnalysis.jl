abstract type AbstractSVDRotation <: AbstractLinearENARotation end
struct SVDRotation <: AbstractSVDRotation end

"""
    SVDRotation()

The default rotation for linear ENA models. Reduces dimensions usings [singular value decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition)
"""
SVDRotation