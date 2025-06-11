#! /bin/env julia


########################################################################
# Computes a packing (i.e., a net) of n points in constant deciderim in
# linear time using hashing, based on Har-Peled and Raichel Net &
# Prune paper.
#
# Implemented by Sariel Har-Peled
########################################################################
#
# 2025-June-7
########################################################################


myBase = "/home/sariel/prog/25/wspd"

push!(LOAD_PATH, myBase * "/src/")
push!(LOAD_PATH, myBase * "/src/cg/")

### To make emacs lsp mode works correctly...
macro ignore(args...) end

@ignore    include("../src/cg/FrechetDist.jl")
@ignore    include("../src/cg/polygon.jl")
@ignore    include("../src/cg/point.jl")
@ignore    include("../src/cg/cg.jl")

using Distributions
using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
using Printf;
using StaticArrays;
using Cairo

include( "graphics.jl" )

@ignore Point{D,T} = FrechetDist.cg.point.Point{D,T}
@ignore Dist =  FrechetDist.cg.point.Dist
@ignore Points =  FrechetDist.cg.polygon.Points

######################################################################
# Misc
######################################################################

@inline function  gid(P::Point{D,T},
                       r::T)::Point{D,Int} where{D,T}
    x = MVector{D,Int}(undef);

    for i ∈ 1:D
        x[i] = floor( Int, P[ i ] / r );
    end

    return  Point{D,Int}( x );
end


######################################################################
# Weighted point set
######################################################################
mutable struct WPoints{D,T}
    orig_PS::Vector{Point{D,T}};    # Original point setup
    PS::Vector{Int};                # The point set -- locations are
                                    # pointers to original points.
    W::Vector{Int};                 # Weights of points
end

@inline function  Base.length( P::WPoints{D,T} ) where {D,T}
    return  length( P.PS );
end


# An empty weighted point set, using the same ground set for defining
# the point set
function  WPoints_empty( G::WPoints{D,T} )  where{D,T}
    return  WPoints( G.orig_PS, Vector{Int}(), Vector{Int}() );
end

function  WPoints( _PS::Vector{Point{D,T}} )  where{D,T}
    PS = [i for i ∈ 1:length(_PS) ]
    W = fill(1, length( PS ) );

    return  WPoints( _PS, PS, W );
end

@inline function  Base.getindex(P::WPoints{D,T}, i::Int) where {D,T}
    return   P.orig_PS[ P.PS[ i ] ];
end


@inline function  weight(P::WPoints{D,T}, i::Int) where {D,T}
    return   P.W[ i ];
end


function   add_weight_to_point( P::WPoints{D,T}, i::Int, w::Int ) where {D,T}
    P.W[ i ] += w;
end

function  oindex( P::WPoints{D,T}, i ) where {D,T}
    return  P.PS[ i ];
end

function  add_point!( P::WPoints{D,T}, oind::Int, w::Int ) where {D,T}
    push!( P.PS, oind );
    push!( P.W , w );
    @assert( length( P.PS ) == length( P.W ) );
    return  length( P.PS );
end


######################################################################
# Linear time algorithm using hashing - using vectors to implement
# open hashing - probably a terrible idea.
######################################################################

"""
GGrid

Generic grid type - every cell stores one value of type VT.
"""
mutable struct GGrid{D,VT}
    ℓ::Float64  # Side length
    cells::Dict{Point{D,Int},VT}
end

function  GGrid_init( ℓ, D, VT )
    return   GGrid{D,VT}( ℓ, Dict{Point{D,Int},VT}() );
end

function  GGrid_store( G::GGrid{D,VT}, trg::Point{D,Int}, v::VT )  where{D,VT}
    G.cells[ trg ] = v;
end

function  GGrid_store( G::GGrid{D,VT}, p::Point{D,T}, v::VT )  where{D,T,VT}
    G.cells[ gid( p, G.ℓ) ] = v;
end

function  GGrid_store( G::GGrid{D,VT}, pnts::WPoints{D,T}, v::VT )  where{D,T,VT}
    for p ∈ pnts
        G.cells[ gid( p, G.ℓ) ] = v;
    end
end


function  GGrid_add_shadow( G::GGrid{D,VT}, p::Point{D,T}, v::VT, nbr )  where{D,T,VT}
    base = gid( p, G.ℓ );
    for  diff ∈ nbr
        trg = base + diff;
        sum = zero( VT );
        if  haskey( G.cells, trg )
            sum = G.cells[ trg ]
        end

        G.cells[ trg ] = sum + v;
    end
end


function  GGrid_add_shadow( G::GGrid{D,VT}, P::WPoints{D,T}, nbr )  where{D,T,VT}
    for i ∈ 1:length( P )
        GGrid_add_shadow( G, P[i], weight( P, i ), nbr );
    end
end



mutable struct WGrid{D,T}
    ℓ::Float64  # Side length
    cells::Dict{Point{D,Int},Vector{Int}}
    P::WPoints{D,T}
end


"""
    store_value!

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

function   WGrid_init( Ground::WPoints{D,T}, ℓ::T )::WGrid{D,T}  where{D,T}
    dict = Dict{Point{D,Int},Vector{Int}}();

    return  WGrid( ℓ, dict, WPoints_empty( Ground ) );
end

function   packing_add_point( G::WGrid{D,T}, loc::Int, rad::T,
                              P::WPoints{D,T} ) where{D,T}
    p =  P[ loc ];
    trg = gid( p, G.ℓ )
    id_min = trg .- 1;
    id_max = trg .+ 1;

    for  cell ∈ CartesianIndex( id_min... ):CartesianIndex( id_max... )
        i = Point{D,Int}( Tuple( cell )... );
        if  ! haskey( G.cells, i )  continue  end
        list = G.cells[ i ];
        for  p_ind ∈ list
            if   Dist( G.P[ p_ind ], p ) < rad
                add_weight_to_point( G.P, p_ind, weight( P, loc ) );
                return
            end
        end
    end

    oloc = oindex( P, loc )
    i = add_point!( G.P, oloc, weight( P, loc ) );
    store_value!( G.cells, trg, i );
end


function  packing( P::WPoints{D,T}, rad::T ) where{D,T}
    G = WGrid_init( P, rad );
    for  i ∈ 1:length( P )
        packing_add_point( G, i, rad, P );
    end

    return  G.P;
end


function  packing( _P::Vector{Point{D,T}}, rad::T ) where{D,T}
    P = WPoints( _P );

    N = packing( P, rad );

    out = Vector{Int}();
    for  i ∈ 1:length(N)
        push!( out, oindex( N, i ) );
    end

    return  out;
end


############################################################################
############################################################################


function   grid_neighberhood( base::Point{D,Int}, dist::T, ℓ::T ) where {D,T}
    Δ = ceil(Int, dist / ℓ ) + 1;

    id_min = base .- Δ
    id_max = base .+ Δ

    #    ℓ = dist / sqrt( D )
    #n = length( P )
    #G = grid_store( P, ℓ, 1:n )
    nbr = Vector{Point{D,Int}}();
    for  _subc ∈ CartesianIndex( id_min... ):CartesianIndex( id_max... )
        subc = Point{D,Int}( Tuple( _subc )... )
        sum::Float64 = 0
        for i ∈ 1:D
            x::Float64 = max( abs( subc[ i ] - base[ i ] ) - 1, 0 );
            sum += x*x;
        end
        d = sqrt(sum) * ℓ;
        if  ( d > dist )
            continue;
        end
        push!( nbr, subc );
    end

    return  nbr;
end


@inline function  check_neighbors( G::WGrid{D,T}, cell, i_pnt,
                                   far, rad, nbr ) where{D,T}
    P = G.P;
    p = P[ i_pnt ];
    for  _subc ∈ nbr
        subc = cell + _subc;
        if  haskey( G.cells, subc )
            subl = G.cells[ subc ]
            for  j ∈ subl
                if  j == i_pnt
                    continue
                end
                if  Dist( p, P[ j ] ) <= rad
                    far[ i_pnt ] = far[ j ]= false
                    return
                end
            end
        end
    end
end

"""
    grid_store

Stor all the points in P[ rng ] in the grid G.
"""

function   grid_store( G, P::WPoints{D,T}, rng ) where{D,T}
    for i ∈ rng
        p = P[ i ]
        trg = gid( p, G.ℓ )
        add_point!( G.P, oindex( P, i ), P.W[ i ] )
        store_value!( G.cells, trg, length( G.P) );
    end
end

function  r_far( P::WPoints{D,T}, rad::T ) where{D,T}
    ℓ = rad / sqrt( D )
    n = length( P )
    G = WGrid_init( P, ℓ )

    nbr = grid_neighberhood( zero(Point{D,Int}), rad, ℓ );

    grid_store( G, P, 1:n )

    #Δ = ceil(Int, sqrt( D ) )

    far = fill( true, n )
    for  (cell, list) ∈ G.cells
        # If there are 2 or more points in the cell, then they are all near...
        if  length( list ) > 1
            for i ∈ list
                far[ i ] = false;
            end
            continue;
        end
        check_neighbors( G, cell, list[ 1 ], far,  rad, nbr );
    end

    return  far;
end

############################################################################
############################################################################


function  nn_dist( P, ind::Int )
    n = length( P );
    @assert( n > 1 );
    p = P[ ind ];
    min_ind = ( ind == 1 ) ? 2 : 1
    min_dist = Dist( p, P[ min_ind ] );
    for  j ∈ 1:n
        if  j == ind  continue  end

        d = Dist( p, P[ j ] )
        if  d < min_dist
            min_dist = d
            min_ind = j
        end
    end

    return  min_ind, min_dist
end

"""
    random_nn_dist

# Returned value

Returns a triple: i, j, dist.

i: the index of a random point in P,
j: The index of the closest point in P to P[i].
dist: Distance between P[i] and P[j].
"""
function  random_nn_dist( P )
    i = rand( 1:length( P ) );
    return  i, nn_dist( P, i )...;
end


function  max_weight( P::WPoints{D,T} ) where {D,T}
    i = argmax( P.W );
    return  i, P.W[ i ]...;
end



function  subset( P::WPoints{D,T}, keep ) where {D,T}
    @assert( length( P ) == length( keep ) );

    O = WPoints_empty( P );

    for  i ∈ 1:length( P.W )
        if   keep[ i ]
            add_point!( O, oindex( P, i ), P.W[ i ] )
        end
    end

    return  O;
end


function   min_ball( P::Vector{Point{D,T}}, i, k ) where {D,T}
    cen = P[ i ];

    arr = [Dist(cen, q) for q ∈ P ];
    d = partialsort!( arr, k )

    return  cen, d;
end

function   min_ball( P::Vector{Point{D,T}}, cen::Point{D,T}, k ) where {D,T}
    arr = [Dist(cen, q) for q ∈ P ];
    d = partialsort!( arr, k )

    return  cen, d;
end

############################################################################
# Grid packing are similar to packing, except that the only constraint
# is that any grid cell has only a single point in the packing. No
# need to check the neighborhood.

function  GridPacking_compute( P::WPoints{D,T}, ℓ::Float64 ) where{D,T}
    G = WGrid_init( P, ℓ );

    for i ∈ 1:length(P)
        trg = gid( P[ i ], G.ℓ )
        if  ! haskey( G.cells, trg )
            add_point!( G.P, oindex( P, i ), weight( P, i ) )
            store_value!( G.cells, trg, length( G.P ) );
            continue;
        end
        list = G.cells[ trg ]
        @assert( length( list ) == 1 );
        ind_p = list[ 1 ]
        G.P.W[ ind_p ] += weight( P, i );
    end
end



############################################################################
############################################################################
"""
    decider_k_ball

# Returned value
  -1 : opt < r
   0 : opt ∈ [r, (1+ε)r]
  +1 : opt > (1+ε)r.
"""
function   decider_k_ball( P::WPoints{D,T}, k::Int, ε::Float64, r::Float64
                         ) where{D,T}
    ℓ = (r/sqrt(D)) * ε / 2.0;
    GS = GGrid_init( ℓ, D, Int );
    nbr = grid_neighberhood( zero(Point{D,Int}), r, ℓ );

    GGrid_add_shadow( GS, P, nbr );
    for  (cell, val) ∈ GS.cells
        if  val < k
            continue
        end
        p = Point{D,T}( cell )
        #println( "PPP : ", p );
        p = ℓ * p;
        #println( "PPPX: ", p );
        return  -1, p;
    end
    return  1, zero( Point{D,T} );
end


function  throw_far_points( P::WPoints{D,T}, r::Float64 )where{D,T}
    far = r_far( P, r );
    near = .!far;
    if ( sum( near ) < length( P ) )
        return   subset( P, near );
    end
    return  P;
end

function smallest_k_ball_binary_search( _P::Vector{Point{D,T}}, k, r, R, ε
                                      )where{D,T}
    VP = WPoints( _P );
    P = packing( VP, ε * r / 4 );

    cen = zero( Point{D,T} );
    iter = 0;

    f_assigned = false;
    while  ( (1+ε)*r < R )  ||  ( iter == 0 )
        iter += 1;
        δ = min( max(  ((R / r ) - 1.0) / 8.0, ε / 4.0 ), 0.5 );
#        if  δ > 1.0
#            δ = 0.5
#        end
        mid = ( r + R ) / 2.0
        res, px = decider_k_ball( P, k, δ, mid );
        if  ( res < 0 )
            f_assigned = true;
            cen = px;
            R = mid;
            continue;
        end
        r = mid;
    end

    if  ! f_assigned
 #       println( "Had to force check..." );
        res, cen = decider_k_ball( P, k, ε, R );
 #       println( "RES= ", res );
        @assert( res < 0 );
    end
    return  min_ball( _P, cen, k );
end

function   smallest_k_ball( _P::Vector{Point{D,T}}, k::Int, ε::Float64 ) where{D,T}
    P = WPoints( _P );

    while  true
        _, _, dist = random_nn_dist( P );

        res, _ = decider_k_ball( P, k, 1.0, dist );

        # Found ball ⟹ Prune: throw away the far points...
        if  ( res < 0 )
            P = throw_far_points( P, dist );
            continue;
        end

        res4, _ = decider_k_ball( P, k, 1.0, 4*dist );

        # Found ball ⟹ optimal solution radius around [r,4r]
        if  ( res4 < 0 )
            return  smallest_k_ball_binary_search( _P, k, dist / 8, 8*dist, ε );
        end

        # Packing (well, net in the original paper)
        N = packing( P, dist );
        mw = max_weight( N )[ 2 ]
        @assert( mw < k );

        P = N;
    end
end


############################################################################
############################################################################

function  draw_two_sets( P, N, filename )
    # Output to file...
    c,cr,_ = cairo_setup( filename, [P], true );

    n = length( P );
    set_source_rgb(cr, 0.8, 0.2, 0.2) # An nice red color
    draw_points( cr, Points( N ), 1.0/(sqrt(n)) );
    Cairo.stroke( cr );

    set_source_rgb(cr, 1.0, 0.0, 0.8);
    draw_points( cr, Points( P ), 1/(16*sqrt(n)) );
    Cairo.stroke( cr );

    finish( c );
    println( "     Created: ", filename );
end


function  test_far_points( P, rad )
    far = r_far(WPoints( Points( P ) ), rad );

    N = Polygon2F();
    for  i ∈ eachindex( P )
        if  ( far[ i ] )
             push!( N, P[ i ] );
        end
    end

    draw_two_sets( P, N, "out/far_points.pdf" );
end


function  test_packing( P, rad )
    net = packing( Points( P ), rad );

    N = Polygon2F();
    for  i ∈ net
        push!( N, P[ i ] );
    end

    draw_two_sets( P, N, "out/net.pdf" );
end



function     test_smallest_ball( P, k )
    cen_b, rad_b = smallest_k_ball( Points( P ), k, 0.5 );
    cen_c, rad_c = smallest_k_ball( Points( P ), k, 0.125 );

    filename = "out/smallest_disk.pdf";
    c,cr,_ = cairo_setup( filename, [P], true );

    set_source_rgba(cr, 0.0, 0.0, 1.0, 0.3) # An nice red color
    Cairo.arc(cr, cen_b[1], cen_b[2], rad_b, 0, 2*pi)
    Cairo.fill(cr)
    Cairo.stroke( cr );

    set_source_rgba(cr, 0.0, 0.0, 0.0, 0.3) # An nice red color
    Cairo.arc(cr, cen_c[1], cen_c[2], rad_c, 0, 2*pi)
    Cairo.fill(cr)
    Cairo.stroke( cr );

    n = length( P );
    set_source_rgb(cr, 0.8, 0.2, 0.2) # An nice red color
    draw_points( cr, Points( P ), 0.1/(sqrt(n)) );
    Cairo.stroke( cr );
    finish( c );

    println( "     Created: ", filename );
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

    println( "\n\n--------------------\n\n" );
    println( "Generating input & Copying..." );
    P = Polygon_random_gaussian( D,Float64, n );

    test_far_points( P, rad );
    test_packing( P, rad );

    k = 20 #div( length( P ), 2);
    k_a = div(k, 3)
    base = npoint( 2.5, 1.5 );
    for i ∈ 1:k_a
        p = Point_random( D, Float64 );
        push!( P, base + p*0.1 );
    end
    test_smallest_ball( P, k )



    return  0;
end
