include("tchakaloff_quad.jl")

maxxlim = 20
nsamp = 1000
nsamp2d = 100

#### Test 1D
d = InducedDistribution(:hermite, 6, Q=500);
d0 = _dists[d.poly]
# Compare PDFs
plot(
    begin
        plot(d0, label=false, c=:black, linewidth=5)
        histogram!(rand(d0, nsamp), bins=50, alpha=0.5, c=:blue, normalize=:pdf, label=false)
        plot!(title="Original Distribution")
        plot!(xlims=(max(-maxxlim,minimum(d0)),min(maxxlim,maximum(d0))))
    end,
    begin
        (d.poly == :laguerre) ? plot(0:0.1:20, x -> pdf(d,x), label=false, c=:black, linewidth=5) : 
            plot(d, label=false, c=:black, linewidth=5)
        histogram!(rand(d, nsamp), bins=50, alpha=0.5, c=:red, normalize=:pdf, label=false)
        plot!(title="Induced Distribution")
        plot!(xlims=(max(-maxxlim,minimum(d0)),min(maxxlim,maximum(d0))))
    end,
    layout=(2,1),
    size=(600,400)
)

## Test 2D
xs = Dict(
    :hermite => range(-8,8,50),
    :legendre => range(-1,1,50),
    :laguerre => range(0,20,50),
    :jacobi2 => range(0,1,50),
    :jacobi4 => range(0,1,50),
    :jacobi09 => range(0,1,50),
)
polys = (:jacobi4,:laguerre)
d = MultivariateInducedDistribution(polys, 5, Q=100);
# Compare PDFs
vissamples = true

plot(
    begin
        heatmap(xs[polys[1]], xs[polys[2]], (x,y) -> pdf(_dists[polys[1]],x)*pdf(_dists[polys[2]],y), title="Original Distribution")
        if vissamples 
            X = reduce(hcat, rand(_dists[d.polys[j]], nsamp2d) for j in 1:2)'
            scatter!(X[1,:], X[2,:], c=:red, label=false)
        else
            plot!()
        end
    end,
    begin
        heatmap(xs[polys[1]], xs[polys[2]], (x,y) -> pdf(d,[x,y]), title="Induced Distribution")
        if vissamples 
            X = rand(d, nsamp2d)
            scatter!(X[1,:], X[2,:], c=:red, label=false)
        else
            plot!()
        end
    end,
    layout=(2,1),
    size=(600,800)
)

nsamp = 100000
plot(
    heatmap(xs[polys[1]], xs[polys[2]], (x,y) -> pdf(_dists[polys[1]],x)*pdf(_dists[polys[2]],y), title="Original Distribution"),
    begin
        X = reduce(hcat, rand(_dists[d.polys[j]], nsamp) for j in 1:2)'
        heatmap(kde((X[1,:], X[2,:])), label=false)
        plot!(xlims=(xs[polys[1]][1],xs[polys[1]][end]),ylims=(xs[polys[2]][1],xs[polys[2]][end]))
        title!("Original Samples")
    end,
    heatmap(xs[polys[1]], xs[polys[2]], (x,y) -> pdf(d,[x,y]), title="Induced Distribution"),
    begin
        X = rand(d, nsamp)
        heatmap(kde((X[1,:], X[2,:])), label=false)
        plot!(xlims=(xs[polys[1]][1],xs[polys[1]][end]),ylims=(xs[polys[2]][1],xs[polys[2]][end]))
        title!("Induced Samples")
    end,
    layout=(2,2),
    size=(800,800)
)

plot(
    begin
        surface(xs[polys[1]], xs[polys[2]], (x,y) -> pdf(_dists[polys[1]],x)*pdf(_dists[polys[2]],y), title="Original Distribution")
        if vissamples 
            X = reduce(hcat, rand(_dists[d.polys[j]], nsamp2d) for j in 1:2)'
            zs = [pdf(_dists[polys[1]],x)*pdf(_dists[polys[2]],y) for (x,y) in zip(X[1,:],X[2,:])]
            scatter!(X[1,:], X[2,:], zs, c=:red, label=false)
        else
            plot!()
        end
    end,
    begin
        surface(xs[polys[1]], xs[polys[2]], (x,y) -> pdf(d,[x,y]), title="Induced Distribution")
        if vissamples 
            X = rand(d, nsamp2d)
            scatter!(X[1,:], X[2,:], pdf.(Ref(d), eachcol(X)), c=:red, label=false)
        else
            plot!()
        end
    end,
    layout=(2,1),
    size=(600,800)
)