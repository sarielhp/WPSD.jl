#! /bin/env julia


push!(LOAD_PATH, pwd()*"/src/")
push!(LOAD_PATH, pwd()*"/src/cg/")

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
#using VirtArray;
using Printf;
using StaticArrays;

mutable struct GridType{D,T}
    P::Vector{Point{D,T}}
    ℓ::Float64  # Side length
    cells::Dict{Point{D,Int},Vector{Int}}

    cp::Tuple{Int, Int};
    cp_dist::T;
    f_regrid::Bool 
end


function  gid(P::Point{D,T}, r::T)::Point{D,Int} where{D,T}
    x = MVector{D,Int}(undef);

    for i ∈ 1:D
        x[i] = floor( Int, P[ i ] / r );
    end

    return  Point{D,Int}( x );
end

"""
    add_value!

    An insertion function for dictionary where a value is an array of values.
"""
function add_value!(dict, key, value)
    if haskey(dict, key)
        push!(dict[key], value)
        #println( "NL: ", length( dict[key] ) );
    else
        dict[key] = [value]
    end
end

function   grid_init( P::Vector{Point{D,T}}, ℓ::T,
                      r::UnitRange{Int} = eachindex( arr ), cp = (-1,-1) ) where{D,T}

    dict = Dict{Point{D,Int},Vector{Int}}();
    sizehint!( dict, min( length( r ), length( P ) ) )
    G = GridType( P, ℓ, dict, cp, ℓ, false );

    for  i ∈ r
        add_value!( G.cells, gid( P[ i ], ℓ ), i );
    end

    return  G;
end

function  closest_pair_add_point( G::GridType{D,T}, loc::Int ) where{D,T}
    p =  G.P[ loc ];
    id = gid( p, G.ℓ )
    v = Point{D,Int}( fill( 1, D ) );
    id_min = sub( id, v );
    id_max = add( id, v );

    #println( "CPAP: ", loc );
    P = G.P;
    #println( "ID: ", id );
    min_ind = -1;
    min_dist = G.cp_dist;
    for  cell ∈ CartesianIndex( id_min... ):CartesianIndex( id_max... )
        i = Point{D,Int}( collect(Tuple( cell ) ) );
        #println( typeof( i ) )
        if  ! haskey( G.cells, i )  continue  end
        list = G.cells[ i ];
        #println( length( list ) );
        if  length( list ) > 4
            G.f_regrid = true;
        end
        #println( "new_dist :", list );
        for  p_ind ∈ list
            new_dist = Dist( P[ p_ind ], p )
            if   new_dist < min_dist
                min_dist = new_dist;
                min_ind = p_ind;
            end
        end
    end

    # Minimum distance had not changed. Nothing much to do...
    if   min_ind < 0
        add_value!( G.cells, gid( P[ loc ], G.ℓ ), loc );
    else
        #println( "min_dist: ", min_dist );
        G.cp = (min_ind, loc );
        G.cp_dist = min_dist;
    end

    if  G.f_regrid
        println( "REGRID!" );
        G = grid_init( G.P, G.cp_dist, 1:loc, G.cp );
    end

    return  G;
end


function  closest_pair( P::Vector{Point{D,T}} ) where{D,T}
    d = Dist( P[ 1 ], P[ 2 ] );
    @assert( d > 0.0 );

    G = grid_init( P, d, 1:2, (1,2) );
    #println( typeof( G ) );
    for  i ∈ 3:length( P )
        G = closest_pair_add_point( G, i );
    end

    return  G;
end

function (@main)(ARGS)
    if  length( ARGS ) != 1
        println( "Usage:\n\t"
                 * "closest_pair.jl [n]\n\n" );
        return  -1;
    end

    D = 2;
    n =  parse(Int64, ARGS[ 1  ] );

    P = Polygon_random( D,Float64, n );

    G = closest_pair( Points( P ) );
    @time G = closest_pair( Points( P ) );
    println( "CP indices:  ", G.cp );
    println( "CP distance: ", G.ℓ );
end
