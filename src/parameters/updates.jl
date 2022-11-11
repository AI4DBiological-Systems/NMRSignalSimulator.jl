
"""
```
updateparameters!(A::CoherencePhase{T})
```

Updates the other fields of A using the contents of `A.κs_β`.
"""
function updateparameters!(A::CoherencePhase{T}, Δc_bar::Vector{Vector{Vector{T}}}) where T <: AbstractFloat
    
    @assert length(Δc_bar) == length(A.κs_β)

    β = NaN

    for i in eachindex(A.κs_β)
        @assert length(Δc_bar[i]) == length(A.cos_β[i]) == length(A.sin_β[i])

        for k in eachindex(A.κs_β[i])
            β = dot(A.κs_β[i], Δc_bar[i][k])
            A.cos_β[i][k] = cos(β)
            A.sin_β[i][k] = sin(β)
        end
    end

    return nothing
end

"""
```
updateparameters!(A::CoherencePhase{T})
```

Updates the other fields of A using the contents of `A.κs_β`.
"""
function updateparameters!(A::CoherenceShift{T}, Δc_bar::Vector{Vector{Vector{T}}}) where T <: AbstractFloat
    
    @assert length(Δc_bar) == length(A.κs_d)

    β = NaN

    for i in eachindex(A.κs_d)
        @assert length(Δc_bar[i]) == length(A.d[i])

        for k in eachindex(A.κs_d[i])
            A.d[i][k] = dot(A.κs_d[i], Δc_bar[i][k])
        end
    end

    return nothing
end