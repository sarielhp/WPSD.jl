#! /bin/env julia


mutable struct GridType{D,T}
    P::Vector{Point{D,T}}
    Dict{Point{D,Int},Vector{Int}}
end 

function add_value!(dict, key, value)
    if haskey(dict, key)
        push!(dict[key], value)
    else
        dict[key] = [value]
    end
end
