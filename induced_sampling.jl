using ClassicalOrthogonalPolynomials
using Distributions
using FastGaussQuadrature
using KernelDensity
using LinearAlgebra
using Plots
using PolyChaos
using ProgressBars
using Random
using Statistics
using StatsPlots

_measures = Dict(
    :legendre => Uniform_11Measure,
    :hermite => GaussMeasure,
    :laguerre => LaguerreMeasure,
    :jacobi2 => () -> Beta01Measure(2,2),
    :jacobi4 => () -> Beta01Measure(4,4),
    :jacobi09 => () -> Beta01Measure(0.9,0.9),
);

_dists = Dict(
    :legendre => Uniform(-1,1),
    :hermite => Normal(0,1),
    :laguerre => Exponential(1.0),
    :jacobi2 => Beta(2,2),
    :jacobi4 => Beta(4,4),
    :jacobi09 => Beta(0.9,0.9),
);

include("index_sets.jl")
include("univariate.jl")
include("multivariate.jl")
include("givens_updowndate.jl")
include("ls_quad.jl")