"""
    Polylogarithmshack

Module containing functions to calculate the polylogarithm and associated functions

"""
module PolylogarithmsHack
#import SpecialFunctions
#import StaticNumbers
#import ForwardDiff
using SpecialFunctions
using StaticNumbers
using ForwardDiff:Dual,value,partials

# not using MPFR for the moment
# using Base.MPFR: ROUNDING_MODE, big_ln2

export polylog, bernoulli, harmonic, stieltjes, dirichlet_beta,polylog_exp,polylog_exp_norm
export Diagnostics
export parse

polylog(x::Real,d::Dual{T,V,N}) where {T,V,N} = Dual{T}(polylog(x,value(d)), polylog(x-1,value(d))/value(d) * partials(d))
#polylog(x::Complex,d::ForwardDiff.Dual{T,V,N}) where {T,V,N} = ForwardDiff.Dual{T}(polylog(x,ForwardDiff.value(d)), polylog(x-1,ForwardDiff.value(d))/ForwardDiff.value(d) * ForwardDiff.partials(d))




# @compat ComplexOrReal{T} = Union{T,Complex{T}}
# s::ComplexOrReal{Float64}
ComplexOrReal{T} = Union{T,Complex{T}}

# Some extra parsing routines for reading Mathematics output, but also, complex numbers
include("utilities.jl")

# Constants
include("constants.jl")

# Series
include("stieltjes.jl")
include("bernoulli_n.jl")
include("harmonic.jl")
# include("gamma_derivatives.jl") # this just generates a table that we include into the code

# Functions
include("beta.jl")
include("bernoulli_poly.jl")
include("polylog.jl")

end
