#! /bin/env julia


########################################################################
# Computes a packing (i.e., a net) of n points in constant dim in
# linear time using hashing, based on Har-Peled and Raichel Net &
# Prune paper.
#
# Implemented by Sariel Har-Peled
########################################################################
#
# 2025-June-7
########################################################################


push!(LOAD_PATH, pwd()*"/src/")
push!(LOAD_PATH, pwd()*"/src/cg/")

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
using Printf;
using StaticArrays;

include( "graphics.jl" )


######################################################################
# Misc
######################################################################

@inline function  gid(P::Point{D,T}, r::T)::Point{D,Int} where{D,T}
    x = MVector{D,Int}(undef);

    for i ∈ 1:D
        x[i] = floor( Int, P[ i ] / r );
    end

    return  Point{D,Int}( x );
end


######################################################################
# Linear time algorithm using hashing - using vectors to implement
# open hashing - probably a terrible idea.
######################################################################


mutable struct GridType{D,T}
    P::Vector{Point{D,T}}
    ℓ::Float64  # Side length
    cells::Dict{Point{D,Int},Vector{Int}}
end


"""
    add_value!

    An insertion function for dictionary where a value is an array of values.
"""
@inline function store_value!(dict::Dict{Point{D,Int},Vector{Int}},
                            key::SVector{D,Int}, value::Int) where {D}
    if haskey(dict, key)
        push!(dict[key], value)
    else
        dict[key] = [value]
    end
end

function   grid_store( P::Vector{Point{D,T}}, ℓ::T,
                      r::UnitRange{Int} = eachindex( P ) ) where{D,T}

    dict = Dict{Point{D,Int},Vector{Int}}();
    sizehint!( dict, min( length( r ), length( P ) ) )
    G = GridType( P, ℓ, dict );

    for  i ∈ r
        trg = gid(P[i], ℓ )
        store_value!( G.cells, trg, i );
    end

    return  G;
end

function   packing_add_point( G::GridType{D,T}, loc::Int, rad::T ) where{D,T}
    p =  G.P[ loc ];
    trg = gid( p, G.ℓ )
    id_min = trg .- 1;
    id_max = trg .+ 1;

    P = G.P;
    for  cell ∈ CartesianIndex( id_min... ):CartesianIndex( id_max... )
        i = Point{D,Int}( Tuple( cell )... );
        #println( typeof( i ) )
        if  ! haskey( G.cells, i )  continue  end
        list = G.cells[ i ];
        if  length( list ) > 4
            G.f_regrid = true;
        end
        for  p_ind ∈ list
            if   Dist( P[ p_ind ], p ) < rad
                return
            end
        end
    end

    store_value!( G.cells, trg, loc );
end


function  packing( P::Vector{Point{D,T}}, rad::T ) where{D,T}
    G = grid_store( P, rad, 1:1);
    for  i ∈ 2:length( P )
        packing_add_point(G, i, rad );
    end

    out = Vector{Int}();
    for  (cell, list) ∈ G.cells
        Base.append!( out, list );
    end

    return  out;
end


############################################################################
############################################################################

@inline function  check_neighbors( P, G, cell, list, far )
    id_min = cell .- Δ
    id_max = cell .+ Δ

    i_pnt = list[ 1 ] ;
    p = P[ i_pnt ];
    for  _subc ∈ CartesianIndex( id_min... ):CartesianIndex( id_max... )
        subc = Point{D,Int}( Tuple( _subc )... )
        if  haskey( G.cells, subc )
            subl = G.cells[ subc ];
            for  j ∈ subl
                if  Dist( p, P[ j ] ) <= rad
                    far[ i_pnt ] = false;
                    return;
                end
            end
        end
    end
end


function  r_far( P::Vector{Point{D,T}}, rad::T ) where{D,T}
    ℓ = rad / sqrt( D )
    n = length( P )
    G = grid_store( P, ℓ, 1:n )

    Δ = ceil(Int, sqrt( D ) )

    out = Vector{Int}()
    far = fill( true, n )
    for  (cell, list) ∈ G.cells
        # If there are 2 or more points in the cell, then they are all near...
        if  length( list ) > 1
            for i ∈ list
                far[ i ] = false;
            end
            continue;
        end
        check_neighbors( P, G, cell, list, far );
    end

    return  far;
end

############################################################################
############################################################################


function  force_compile()
    return  0
end


function (@main)(ARGS)
    if  length( ARGS ) != 2
        println( "Usage:\n\t"
                 * "packing.jl [n] [rad]\n\n" );
        return  -1;
    end

    D = 2;
    n =  parse(Int64, replace( ARGS[ 1  ], "," => "" ) );
    rad =  parse(Float64, replace( ARGS[ 2  ], "," => "" ) );

    println( "n: ", n );
    println( "rad: ", rad );
    force_compile();

    println( "\n\n--------------------\n\n" );
    println( "Generating input & Copying..." );
    P = Polygon_random_gaussian( D,Float64, n );

    net = r_far(Points( P ), rad );

    # Output to file...
    c,cr,bb = cairo_setup( "out.pdf", [P], true );

    N = Polygon2F();
    for  i ∈ net
        #println( "i: ", i );
        push!( N, P[ i ] );
    end

    println( "LEN: ", length( P ) );
    println( "N LEN: ", length( N ) );
    set_source_rgb(cr, 0.8, 0.2, 0.2) # An nice red color
    draw_points( cr, Points( N ), 0.3/(sqrt(n)) );
    Cairo.stroke( cr );


    set_source_rgb(cr, 1.0, 0.0, 0.8);
    draw_points( cr, Points( P ), 1/(16*sqrt(n)) );
    Cairo.stroke( cr );

    finish( c );

    return  0;
end
