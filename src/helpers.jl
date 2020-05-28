# TODO verify
function help_deflating_svd(networkModel::DataFrame, unitModel::DataFrame, controlModel::Union{Nothing,DataFrame}=nothing)
    # mean center X
    rawCounts = unitModel[!, networkModel[!, :relationship]]
    countMeans = collect(mean(rawCounts[!, i]) for i in 1:size(rawCounts)[2])
    X = Matrix{Float64}(DataFrame(rawCounts .- transpose(countMeans)))

    # for each column in the controls, and for each column in the unit data,
    # deflate the unit data by subtracting out of it the projection of itself
    # onto that control column
    if !isnothing(controlModel)

        # mean center Y
        controlMeans = collect(mean(controlModel[!, i]) for i in 1:size(controlModel)[2])
        Y = Matrix{Float64}(DataFrame(controlModel .- transpose(controlMeans)))
        for i in 1:size(X)[2]
            xcol = X[:, i]
            for j in 1:size(Y)[2]
                v = Y[:, j]
                scalar = dot(xcol, v) / dot(v, v)
                xcol -= scalar * v
            end

            X[:, i] = xcol
        end
    end

    # then, once we've deflated or not, we run an SVD on the data
    Z = transpose(X) * X
    pcaModel = projection(fit(PCA, Matrix{Float64}(Z), pratio=1.0))
    return pcaModel
end

# # see https://www.pnas.org/content/113/51/14662
# # If you only have one column in the controlModel, then you really should use an ortho svd
# # This uses the AC-PCA technique to find a dimensions of the networkModel space that maximize variance,
# # penalized for colinearity with the controlModel, and attempting to optimize orthogonality with the
# # first column of the controlModel. If no controlModel is given, then it just runs an svd.
# function help_ac_svd(networkModel::DataFrame, unitModel::DataFrame, controlModel::Union{Nothing,DataFrame}=nothing, threshold::Real=0.15)
#     rawCounts = unitModel[!, [networkRow[:relationship] for networkRow in eachrow(networkModel)]]
#     countMeans = collect(mean(rawCounts[!, i]) for i in 1:size(rawCounts)[2])
#     X = Matrix{Float64}(DataFrame(rawCounts .- transpose(countMeans)))
#     Z = transpose(X) * X
#     if !isnothing(controlModel)
#         x_axis = controlModel[:, 1]
#         norm_x = sqrt(sum(x_axis .* x_axis))
#         controlMeans = collect(mean(controlModel[!, i]) for i in 1:size(controlModel)[2])
#         Y = Matrix{Float64}(DataFrame(controlModel .- transpose(controlMeans)))
#         K = Y * transpose(Y)
#         hi_lambda = 1
#         lo_lambda = 1
#         best_angle = 1.0
#         best_model = Real[]
#         for i in 1:10
#             Z = transpose(X) * X - hi_lambda * transpose(X) * K * X
#             pcaModel = projection(fit(PCA, Matrix{Float64}(Z), mean=0, pratio=1.0))
#             y_axis_p = map(eachrow(unitModel)) do unitRow
#                 return sum(
#                     pcaModel[i, 1] * unitRow[networkRow[:relationship]]
#                     for (i, networkRow) in enumerate(eachrow(networkModel))
#                 )
#             end

#             norm_y = sqrt(sum(y_axis_p .* y_axis_p))
#             cos_angle = sum(x_axis .* y_axis_p) / norm_x / norm_y
#             if abs(cos_angle) < threshold # close enough, return
#                 return pcaModel
#             elseif abs(cos_angle) < best_angle # we're getting better, keep going
#                 lo_lambda = hi_lambda
#                 hi_lambda *= 2
#                 best_angle = abs(cos_angle)
#                 best_model = pcaModel
#             else # we overshot it, move on to stage two
#                 break
#             end
#         end

#         for i in 1:10
#             lambda = (lo_lambda + hi_lambda) / 2
#             Z = transpose(X) * X - lambda * transpose(X) * K * X
#             pcaModel = projection(fit(PCA, Matrix{Float64}(Z), mean=0, pratio=1.0))
#             y_axis_p = map(eachrow(unitModel)) do unitRow
#                 return sum(
#                     pcaModel[i, 1] * unitRow[networkRow[:relationship]]
#                     for (i, networkRow) in enumerate(eachrow(networkModel))
#                 )
#             end

#             norm_y = sqrt(sum(y_axis_p .* y_axis_p))
#             cos_angle = sum(x_axis .* y_axis_p) / norm_x / norm_y
#             if abs(cos_angle) < threshold # close enough, return
#                 return pcaModel
#             elseif abs(cos_angle) < best_angle # try bigger lambda
#                 lo_lambda = lambda
#                 best_angle = abs(cos_angle)
#                 best_model = pcaModel
#             else # try smaller lambda
#                 hi_lambda = lambda
#             end
#         end
        
#         return best_model # we gave it our best shot
#     else
#         pcaModel = projection(fit(PCA, Matrix{Float64}(Z), mean=0, pratio=1.0))
#         return pcaModel
#     end
# end