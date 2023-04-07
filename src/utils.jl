function derivedAnyCode!(data, newCol, oldCols...)
    data[!, newCol] = ones(nrow(data))
    for col in oldCols
        data[!, newCol] = data[!, newCol] .* (1 .- data[!, col])
    end

    data[!, newCol] = 1 .- data[!, newCol]
end

function derivedAllCode!(data, newCol, oldCols...)
    data[!, newCol] = ones(nrow(data))
    for col in oldCols
        data[!, newCol] = data[!, newCol] .* data[!, col]
    end
end