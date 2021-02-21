function help_deflating_svd(networkModel::DataFrame, unitModel::DataFrame, controlModel::Union{Nothing,DataFrame}=nothing)
    X = Matrix{Float64}(unitModel[!, networkModel[!, :relationship]])
    for i in 1:size(X)[2]
        xcol = X[:, i]
        xcol = xcol .- mean(xcol) # mean center
        X[:, i] = xcol
    end

    if !isnothing(controlModel)
        C = Matrix{Float64}(controlModel)
        for i in 1:size(C)[2]
            ccol = C[:, i]
            ccol = ccol .- mean(ccol) # mean center
            C[:, i] = ccol
            for j in 1:size(X)[2] # deflate
                xcol = X[:, j]
                scalar = dot(xcol, ccol) / dot(ccol, ccol)
                xcol -= scalar * ccol
                X[:, j] = xcol
            end
        end
    end

    # then, once we've deflated or not, we run an SVD on the data
    pcaModel = fit(PCA, X', pratio=1.0)
    return pcaModel
end

function help_plot_ci(p, xs, ys, color, shape, label)
    if length(xs) > 0
        x = mean(xs)
        y = mean(ys)
        Plots.plot!(p, [x], [y], 
            label=label,
            seriestype=:scatter,
            markersize=4,
            markershape=shape,
            markercolor=color,
            markerstrokecolor=color)
    end

    if length(xs) > 1
        ci_x = collect(confint(OneSampleTTest(xs)))
        ci_y = collect(confint(OneSampleTTest(ys)))    
        Plots.plot!(p, [ci_x[1], ci_x[2]], [ci_y[1], ci_y[1]], 
            label=nothing,
            seriestype=:line,
            linewidth=1,
            linecolor=color)

        Plots.plot!(p, [ci_x[1], ci_x[2]], [ci_y[2], ci_y[2]], 
            label=nothing,
            seriestype=:line,
            linewidth=1,
            linecolor=color)

        Plots.plot!(p, [ci_x[1], ci_x[1]], [ci_y[1], ci_y[2]], 
            label=nothing,
            seriestype=:line,
            linewidth=1,
            linecolor=color)

        Plots.plot!(p, [ci_x[2], ci_x[2]], [ci_y[1], ci_y[2]], 
            label=nothing,
            seriestype=:line,
            linewidth=1,
            linecolor=color)
    end
end

function help_nonlinear_gradient(lo, mid, hi; grains=100, curve=1.5)
    return vcat(
        [weighted_color_mean((100-i)^curve/grains^curve, lo, mid) for i in 1:grains],
        [weighted_color_mean(1-i^curve/grains^curve, mid, hi) for i in 1:grains]
    )
end