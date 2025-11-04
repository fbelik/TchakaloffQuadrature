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
d = MultivariateInducedDistribution((:hermite,:laguerre), N, index_set=multi_index_set(2, N, HC), Q=500);

q = least_squares_quad(d, M);
q0 = least_squares_quad(d, M, induced=false);
plot(q, zcolor=q.weights .>= 0)
plot(q0, zcolor=q0.weights .>= 0)

plot(q)
plot!(q0)

succ, plt = trial_comparison(d, M, trials, tol=1e-4); plt

# Continuous insertion/removal
tol=0;
d = InducedDistribution(:hermite,10);
d = MultivariateInducedDistribution((:legendre,:hermite),6);

plt1, plt2, ws, cs, Ms = continuous_insertion_trials(d, tol=tol, MMax=1000);
plt3, plt4, w0s, c0s, M0s = continuous_insertion_trials(d, induced=false, MMax=1000,tol=tol);
plt5, plt6, wrs, crs, Mrs = continuous_removal_trials(d, tol=tol, MMax=1000);
plt7, plt8, w0rs, c0rs, M0rs = continuous_removal_trials(d, induced=false, MMax=1000, tol=tol);

maxxlim = 1000
ylims=(-1e-2,0)
ylims2=(0,5)
mean([w[end] for w in ws])
mean([w[end] for w in wrs])

plot(
    plot!(plt1,title="Induced minimum weight",xlims=(Ms[1][1],maxxlim), ylims=ylims),
    plot!(plt2,title="Induced concentration",xlims=(Ms[1][1],maxxlim), ylims=ylims2),
    plot!(plt3,title="Standard minimum weight",xlims=(Ms[1][1],maxxlim), ylims=ylims),
    plot!(plt4,title="Standard concentration",xlims=(Ms[1][1],maxxlim), ylims=ylims2),
    layout=(2,2),
    size=(600,600),
    plot_title="Continuous Insertion"
)

plot(
    plot!(plt5,title="Induced minimum weight",xlims=(Ms[1][1],maxxlim), ylims=ylims),
    plot!(plt6,title="Induced concentration",xlims=(Ms[1][1],maxxlim), ylims=ylims2),
    plot!(plt7,title="Standard minimum weight",xlims=(Ms[1][1],maxxlim), ylims=ylims),
    plot!(plt8,title="Standard concentration",xlims=(Ms[1][1],maxxlim), ylims=ylims2),
    layout=(2,2),
    size=(600,600),
    plot_title="Continuous Removal"
)

d = InducedDistribution(:hermite, 10);

animnum = 10000
eachframe = floor(Int, animnum / 250)
q = least_squares_quad(d, induced=true); anim = @animate for _ in 1:animnum 
    add_node!(q)
    plot(
        scatter(q.nodes, q.weights, label=false, zcolor=(q.weights .> 0), title="LS Weights", ylims=(min(0,2*minimum(q.weights)),2*maximum(q.weights))),
        scatter(q.nodes, monte_carlo_weights(q), label=false, c=:orange, title="CS Shifted MC Weights", ylims=(0,2*maximum(monte_carlo_weights(q)))),
        scatter(q.nodes, monte_carlo_weights(q, false), label=false, c=:green, title="Original MC Weights", ylims=(0,2*maximum(monte_carlo_weights(q,false)))),
        scatter(q.nodes, q.CS_rescale, label=false, c=:red, title="CS Rescaling Factor", ylims=(0,2*maximum(q.CS_rescale))),
        layout=(2,2),
        size=(600,600),
        plot_title="M=$(q.Mtot), Positive=$(is_positive(q))"
    )
end every eachframe; gif(anim)

q = least_squares_quad(d, induced=true); anim = @animate for _ in 1:animnum 
    remove_add_node!(q)
    plot(
        scatter(q.nodes, q.weights, label=false, zcolor=(q.weights .> 0), title="LS Weights", ylims=(min(0,2*minimum(q.weights)),2*maximum(q.weights))),
        scatter(q.nodes, monte_carlo_weights(q), label=false, c=:orange, title="CS Shifted MC Weights", ylims=(0,2*maximum(monte_carlo_weights(q)))),
        scatter(q.nodes, monte_carlo_weights(q, false), label=false, c=:green, title="Original MC Weights", ylims=(0,2*maximum(monte_carlo_weights(q,false)))),
        scatter(q.nodes, q.CS_rescale, label=false, c=:red, title="CS Rescaling Factor", ylims=(0,2*maximum(q.CS_rescale))),
        layout=(2,2),
        size=(600,600),
        plot_title="M=$(q.Mtot), Positive=$(is_positive(q))"
    )
end every eachframe; gif(anim)
