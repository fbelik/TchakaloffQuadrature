include("induced_sampling.jl")

# Quadrature rule for integation
N = 10; M=200;trials=1000;
d = InducedDistribution(:legendre, N);
d = InducedDistribution(:hermite, N);
d = InducedDistribution(:jacobi09, N);

q = least_squares_quad(d, M);
q0 = least_squares_quad(d, M, induced=false);
plot(q)
plot!(q0)

succ, plt = trial_comparison(d, M, trials, tol=1e-4); plt

# Quadrature rule for multivariate integation
N = 15; M=500;trials=1000;
d = MultivariateInducedDistribution((:hermite,:laguerre), N, index_set=multi_index_set(2, N, :hc), Q=1000);

q = least_squares_quad(d, M);
q0 = least_squares_quad(d, M, induced=false);
plot(q, zcolor=q.weights .>= 0)
plot(q0, zcolor=q0.weights .>= 0)

plot(q)
plot!(q0)

succ, plt = trial_comparison(d, M, trials, tol=1e-4); plt

# Continuous insertion/removal
tol=1e-6;
d = InducedDistribution(:hermite,15);
d = MultivariateInducedDistribution((:legendre,:hermite),6);

plt1, plt2, ws, cs, Ms = continuous_insertion_trials(d, tol=tol, MMax=1000);
plt3, plt4, w0s, c0s, M0s = continuous_insertion_trials(d, induced=false, MMax=1000,tol=tol);
plt5, plt6, wrs, crs, Mrs = continuous_removal_trials(d, tol=tol, MMax=1000);
plt7, plt8, w0rs, c0rs, M0rs = continuous_removal_trials(d, induced=false, MMax=1000, tol=tol);

maxxlim = 1000

plot(
    plot!(plt1,title="Induced minimum weight",xlims=(Ms[1][1],maxxlim)),
    plot!(plt2,title="Induced concentration",xlims=(Ms[1][1],maxxlim)),
    plot!(plt3,title="Standard minimum weight",xlims=(Ms[1][1],maxxlim)),
    plot!(plt4,title="Standard concentration",xlims=(Ms[1][1],maxxlim)),
    layout=(2,2),
    size=(600,600),
    plot_title="Continuous Insertion"
)

plot(
    plot!(plt5,title="Induced minimum weight",xlims=(Ms[1][1],maxxlim)),
    plot!(plt6,title="Induced concentration",xlims=(Ms[1][1],maxxlim)),
    plot!(plt7,title="Standard minimum weight",xlims=(Ms[1][1],maxxlim)),
    plot!(plt8,title="Standard concentration",xlims=(Ms[1][1],maxxlim)),
    layout=(2,2),
    size=(600,600),
    plot_title="Continuous Removal"
)