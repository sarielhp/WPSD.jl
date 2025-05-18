module  WSPD


push!(LOAD_PATH, pwd()*"/cg/")
push!(LOAD_PATH, pwd() )

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
using VirtArray;
using Cairo, Colors

#include("./VirtArray.jl")
#include("./BBT.jl")

using BBT
using DataStructures

mutable struct WSPair{D,T}
    left::BBTNode
    right::BBTNode

    dist::T # Distance between the boxes.
    max_sides_diam::T
    diam_ub::T  # Upper bound on the maximum diameter of any pair of
                # points in this pair.
end

function  WPDPair_init( left::BBTNode{D,T}, right::BBTNode{D,T} )  where {D,T}
    @assert( !isnothing( left ) );
    @assert( !isnothing( right ) );

    dist = BBox_dist( left.bb, right.bb );
    max_sides_diam = max( BBox_diam( left.bb ), BBox_diam( right.bb ) );
    diam_ub = BBox_max_dist( left.bb, right.bb );

    return  WSPair{D,T}( left, right, dist, max_sides_diam, diam_ub );
end

const MAX_SEP = 100000000.0
function  separation( pair::WSPair{D,T} ) where {D,T}
    if  pair.dist <= 0.0
        return  MAX_SEP   # Or should I return +∞
    end
    return  pair.max_sides_diam / pair.dist;
end


#mutable struct  WSPDOrder{T} <: Base.Order.Ordering
#
#end

mutable struct  WSPDOrder{D,T} <: Base.Order.Ordering
    dummy::Int64
end

#function  lt(o::EventsOrder, a::Int64, b::Int64 )
#    return  isless( o.values[ a ], o.values[ b ] );
#end

import Base.Order.lt
function  lt( o::WSPDOrder{D,T}, a::WSPair{D,T}, b::WSPair{D,T} ) where {D,T}
    o.dummy += 1;
    return  a.diam_ub > b.diam_ub;
end

PairInt = Tuple{Int,Int};
mutable struct PD{D,T}
    pnts::Polygon{D,T};
    finals::Vector{WSPair{D,T}};
    curr_active::Vector{WSPair{D,T}};
    heap::BinaryHeap{WSPair{D,T}};
    tree::BBTree{D,T}
    sep::T;  # Desired quality of separation
    hash_pairs::Dict{PairInt, Bool};
end


@inline function  get_id( u::BBTNode, v::BBTNode )::PairInt
    x, y = u.id, v.id;

    x = min( u.id, v.id );
    y = max( u.id, v.id );

    return  (x,y)
end


@inline function  has_pair( W::PD{D,T}, pair::PairInt ) where {D,T}
    return  haskey( W.hash_pairs, pair );
end


function  push_pair( W::PD{D,T}, u::BBTNode{D,T}, v::BBTNode{D,T}
                   ) where {D,T}
    id = get_id( u, v );
    if  has_pair( W, id )
        #println( "REJECTED:", id );
        return
    end;

    #println( id );
    p_a = WPDPair_init( u, v );
    push!( W.heap, p_a );
    W.hash_pairs[ id ] = true;

    return  p_a;
end

"""
    WSPD_top_refine

    Takes the top of the active pairs, and refine it.
"""
function  WSPD_top_refine( W::PD{D,T} ) where {D,T}
    if ( isempty( W.heap ) ) return  end

    top = pop!( W.heap );
    BBTree_refine_node( W.tree, top.left );
    BBTree_refine_node( W.tree, top.right );

    l_diam = BBox_diam( top.left.bb );
    r_diam = BBox_diam( top.right.bb );

    # Is this pair already between points?
    if  ( ( r_diam == 0.0 )  &&  ( l_diam == 0.0 ) )
        push!( W.finals, top );
        return;
    end

    f_refine_left = (l_diam >= r_diam );

    if  f_refine_left
        push_pair( W, top.left.left, top.right );
        push_pair( W, top.left.right, top.right );
    else
        push_pair( W, top.left, top.right.left );
        push_pair( W, top.left, top.right.right );
    end
end


"""
    WSPD_top_delete

    Takes the top of the active pairs, and deletes it.
"""
function  WSPD_top_delete( W::PD{D,T} )  where{D,T}
    if ( isempty( active ) ) return  end

    top = pop!( W.heap );
end


"""
    WSPD_top_finalize

Takes the top of the active pairs, and move it to the generated list of pairs.
"""
function  WSPD_top_finalize( W::PD{D,T} )  where  {D,T}
    if ( isempty( W.heap ) ) return  end
    top = pop!( W.heap );
    push!( W.finals, top );
end


function  WSPD_init( _pnts::Polygon{D,T}, _sep::T ) where {D,T}
    tree = BBTree_init( _pnts );

    PairT = WSPair{D,T}
    #order = Base.Order.ForwardOrdering{PairT}( compare_pairs );
    order =  WSPDOrder{D,T}( 0 );

    curr_active = Vector{PairT}();
    heap = BinaryHeap{PairT}( order, curr_active );
    W = PD( _pnts, Vector{WSPair{D,T}}(),
        curr_active, heap, tree, _sep,
        Dict{PairInt, Bool}() );

    pair = WSPair{D,T}( tree.root, tree.root, 0.0, BBox_diam( tree.root.bb ),
                        BBox_diam( tree.root.bb ) );
    push!( W.heap, pair );

    return  W;
end


function   WSPD_expand( W::PD{D,T} ) where {D,T}
    while  ( !isempty( W.heap ) )
        top = first( W.heap );
        #println( "    Δ: ", top.max_diam, "    [", separation( top ), "]:",
        #    length( W.finals ) );
        if  ( separation( top ) > W.sep )
            WSPD_top_refine( W );
            continue;
        end
        #println( top.left.r, " × ", top.right.r );
        WSPD_top_finalize( W );
    end
end

function   WSPD_get_top( W::PD{D,T} ) where {D,T}
    @assert( ! isempty( W.heap ) );

    return  first( W.heap );
end

function   WSPD_top_diam_ub( W::PD{D,T} ) where {D,T}
    @assert( ! isempty( W.heap ) );
    if isempty( W.heap )
        return  zero( T );
    end

    return  first( W.heap ).diam_ub;
end

function  WSPD_get_reps( W::PD{D,T}, pair::WSPair{D,T} ) where {D,T}
    @assert( ! isempty( pair.left.r ) )
    @assert( ! isempty( pair.right.r ) )
    return  W.tree.PS[ first( pair.left.r ) ], W.tree.PS[ first( pair.right.r ) ]
end

function  WSPD_reps_orig_indexes( W::PD{D,T}, pair::WSPair{D,T} ) where {D,T}
    @assert( ! isempty( pair.left.r ) )
    @assert( ! isempty( pair.right.r ) )

    l_i = first( pair.left.r )
    r_i = first( pair.right.r )

    #println( typeof( W.tree.PS ) );

    #println( orig_index( W.tree.PS, l_i ) );
    return  orig_index( W.tree.PS, l_i ),orig_index( W.tree.PS, r_i )
end


### The exports...

export WPD
export WPDPair

export WSPD_top_delete, WSPD_top_finalize, WSPD_top_delete
export WSPD_expand, WSPD_init;

export WSPD_top_diam_ub
export WSPD_get_top
export WSPD_top_refine

export WSPD_get_reps
export WSPD_reps_orig_indexes 
end
