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

@inline function  tgid(P::Point{D,T}, r::T)::Point{D,Int} where{D,T}
    x = MVector{D,Int}(undef);

    for i ∈ 1:D
        x[i] = floor( Int, P[ i ] / r );
    end

    return  Point{D,Int}( x );
end
@inline function  tgid_trg(P::Point{D,T}, r::T, trg::MVector{D,Int} ) where{D,T}
    for i ∈ 1:D
        trg[i] = floor( Int, P[ i ] / r );
    end
end


######################################################################
# Linear time algorithm using hashing - using vectors to implement
# open hashing - probably a terrible idea.
######################################################################


mutable struct GridType{D,T}
    P::Vector{Point{D,T}}
    ℓ::Float64  # Side length
    cells::Dict{Point{D,Int},Vector{Int}}

    cp::Tuple{Int, Int};
    cp_dist::T;
    f_regrid::Bool
end



"""
    add_value!

    An insertion function for dictionary where a value is an array of values.
"""
@inline function add_value!(dict, key, value)
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
    add_value!( G.cells, gid( P[ loc ], G.ℓ ), loc );
    if   min_ind > 0
        #G.f_regrid = true;
        #println( "min_dist: ", min_dist );
        G.cp = (min_ind, loc );
        G.cp_dist = min_dist;
    end

    if  G.f_regrid
        #println( "REGRID!" );
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


############################################################################
# closest pair using open hashing implemented using an extra counter
############################################################################


mutable struct OHGridType{D,T,DPONE}
    P::Vector{Point{D,T}}
    ℓ::Float64  # Side length
    cells::Dict{Point{DPONE,Int},Int}

    cp::Tuple{Int, Int};
    cp_dist::T;
    f_regrid::Bool
end

function  OHGridType( _P::Vector{Point{D,T}}, _ℓ::Float64,
                      _cells::Dict{Point{DPONE,Int},Int},
                     _cp::Tuple{Int, Int} ) where{D,T,DPONE}
    @assert( (D+1) == DPONE );
    return OHGridType( _P, _ℓ, _cells, _cp, _ℓ, false );
end


"""
    OH_add_value!

    An insertion function for dictionary where a value is an array of values.
"""
@inline function OH_add_value!(dict, _key::Point{D,Int}, value, count = 1) where {D}
    key = Point{D+1,Int}( _key..., count );
    while  haskey( dict, key )
        count += 1;
        key = Point{D+1,Int}( _key..., count );
    end

    dict[key] = value
end


function   OH_extract_list( G::OHGridType{D,T,DPONE}, p::Point{D,Int},
                            out::Vector{Int} ) where{D,T,DPONE}
    empty!( out );
    count = 1;
    key = Point{D+1,Int}( p..., count );
    while  true
        push!( out, G.cells[ key ] );

        count += 1;
        key = Point{D+1,Int}( p..., count );

        if  ! haskey( G.cells, key )  break  end
    end
end


function   OH_grid_init( P::Vector{Point{D,T}}, ℓ::T,
                         r::UnitRange{Int} = eachindex( arr ), cp = (-1,-1) ) where{D,T}

    dict = Dict{Point{D+1,Int},Int}();
    sizehint!( dict, min( length( r ), length( P ) ) )
    G = OHGridType( P, ℓ, dict, cp );

    for  i ∈ r
        OH_add_value!( G.cells, gid( P[ i ], ℓ ), i );
    end

    return  G;
end

function  OH_cp_add_point( G::OHGridType{D,T}, loc::Int, list::Vector{Int} ) where{D,T}
    p =  G.P[ loc ];
    id = gid( p, G.ℓ )
    cid = CartesianIndex( id... );
    v = Point{D,Int}( fill( 1, D ) );
    id_min = sub( id, v );
    id_max = add( id, v );

    P = G.P;

    min_ind = -1;
    min_dist = G.cp_dist;
    id_load::Int = 0;
    for  cell ∈ CartesianIndex( id_min... ):CartesianIndex( id_max... )
        i = Point{D+1,Int}( Tuple( cell )..., 1 );
        if  ! haskey( G.cells, i )  continue  end

        OH_extract_list( G, Point{D,Int}( collect(Tuple( cell ) ) ), list );

        ll = length( list )
        if  ( cell == cid )
            id_load = ll
        end
        if  ll > 4
            G.f_regrid = true;
        end
        for  p_ind ∈ list
            new_dist = Dist( P[ p_ind ], p )
            if   new_dist < min_dist
                min_dist = new_dist;
                min_ind = p_ind;
#                G.f_regrid = true;
            end
        end
    end

    OH_add_value!( G.cells, gid( P[ loc ], G.ℓ ), loc, id_load+1 );
    if   min_ind > 0
        G.cp = (min_ind, loc );
        G.cp_dist = min_dist;
    end

    if  G.cp_dist < G.ℓ  &&  G.f_regrid
#    if  G.f_regrid
        #println( "REGRID!" );
        G = OH_grid_init( G.P, G.cp_dist, 1:loc, G.cp );
    end

    return  G;
end


function  OH_closest_pair( P::Vector{Point{D,T}} ) where{D,T}
    d = Dist( P[ 1 ], P[ 2 ] );
    @assert( d > 0.0 );

    list = Vector{Int}();

    G = OH_grid_init( P, d, 1:2, (1,2) );
    #println( typeof( G ) );
    for  i ∈ 3:length( P )
        G = OH_cp_add_point( G, i, list );
    end

    return  G;
end




#############################################################################
############################################################################
# closest pair using open hashing implemented using preallocated array for
#  the conflict lists
############################################################################

CL_MAX_SIZE::Int64 = 5

mutable struct CLGridType{D,T}
    P::Vector{Point{D,T}}
    ℓ::Float64  # Side length
    cells::Dict{Point{D,Int},Int}

    lists::Vector{Int};
    lists_pos::Int64

    cp::Tuple{Int, Int};
    cp_dist::T;
    f_regrid::Bool
    one::Point{D,Int}
    trg::MVector{D,Int}
end

function  CLGridType( _P::Vector{Point{D,T}}, _ℓ::Float64,
    _cp::Tuple{Int, Int},
    _cells::Dict{Point{D,Int},Int},
    _lists::Vector{Int}
) where{D,T}
    return CLGridType( _P, _ℓ, _cells,
        _lists,
        1,
        _cp, _ℓ, false,
        Point{D,Int}( fill( 1, D ) ),
        MVector{D,Int}( fill( 0, D ) ) );
end


"""
    CL_add_value!

    An insertion function for dictionary where a value is an array of values.
"""
@inline function CL_add_value!( G::CLGridType{D,T}, key::MVector{D,Int},
    value::Int ) where {D,T}
    if ! haskey( G.cells, key )
        loc = G.lists_pos
        G.lists_pos += CL_MAX_SIZE;
        G.cells[ key ] = loc;
    else
        loc = G.cells[ key ];
    end
    start = loc;
    while  ( G.lists[ loc ] != 0 )
        loc += 1
    end
    G.lists[ loc ] = value;
    if  ( loc - start ) > ( CL_MAX_SIZE - 1 )
        G.f_regrid = true;
    end
end


function   CL_grid_init( P::Vector{Point{D,T}}, ℓ::T,
    _cells::Dict{Point{D,Int},Int},
    _lists::Vector{Int},
        r::UnitRange{Int} = eachindex( arr ),
    cp = (-1,-1)
) where{D,T}
    #=
    for  ( key, val) ∈ _cells
        for  i ∈ 0:CL_MAX_SIZE
            _lists[ val + i ] = 0;
        end
    end
    =#
    #@time
    fill!( _lists, zero( Int ) );
    empty!( _cells );
    @assert( ℓ > 0.0 );

    G = CLGridType( P, ℓ, cp, _cells, _lists );

    G.trg = MVector{D,Int}( undef )
    for  i ∈ r
        tgid_trg( P[ i ], ℓ, G.trg )
        CL_add_value!( G, G.trg, i );
    end

    return  G;
end

@inline function  CL_cp_add_point( G::CLGridType{D,T}, loc::Int,
) where{D,T}
    p =  G.P[ loc ];
    id = gid( p, G.ℓ )
    id_min = sub( id, G.one );
    id_max = add( id, G.one );

    P = G.P;

    min_ind = -1;
    min_dist = G.cp_dist;

    for  cell ∈ CartesianIndex( id_min... ):CartesianIndex( id_max... )
        i = Point{D,Int}( Tuple( cell )... );
        if  ! haskey( G.cells, i )  continue  end
        cloc = G.cells[ i ] ;

        while  G.lists[ cloc ] != 0
            p_ind = G.lists[ cloc ] ;
            new_dist = Dist( P[ p_ind ], p )
            if   new_dist < min_dist
                min_dist = new_dist;
                min_ind = p_ind;
            end
            cloc += 1
        end
    end

    tgid_trg( P[ loc ], G.ℓ, G.trg )
    CL_add_value!( G, G.trg, loc );
    if   min_ind > 0
        G.cp = (min_ind, loc );
        G.cp_dist = min_dist;
    end

    return  G;
end


function  CL_closest_pair( P::Vector{Point{D,T}} ) where{D,T}
    d = Dist( P[ 1 ], P[ 2 ] );
    @assert( d > 0.0 );

    cells = Dict{Point{D,Int},Int}();
    sizehint!( cells, 4*length( P ) );
#    println( "len P: ", length( P ) );
#    println( length( cells.slots ) );
#    println( length( cells.keys ) );

#    exit( -1 );
    
    lists = zeros(Int, (2+length( P )) * CL_MAX_SIZE )

    
    G = CL_grid_init( P, d, cells, lists, 1:2, (1,2) );
    for  i ∈ 3:length( P )
        CL_cp_add_point( G, i );

        if  G.cp_dist < G.ℓ  &&  G.f_regrid
            G = CL_grid_init( G.P, G.cp_dist, G.cells, G.lists, 1:i, G.cp );
        end
    end

    return  G;
end





############################################################################
# Closest pair in the plane using divide-and-conquer in O( n log n ) time.
############################################################################


mutable struct  Solution
    d::Float64
    i::Int64
    j::Int64
end
function  Sol_set( sol::Solution, _d, _i, _j )
    sol.d, sol.i, sol.j = _d, _i, _j;
end

Points2F = Vector{Point2F};

struct  PointSetXY
    P::Points2F;
    ord_x::Vector{Int64};
    ord_y::Vector{Int64};
end

len( PS::PointSetXY ) = length( PS.ord_x );    # number of points in PS
#d( p, q ) =  norm( ( p[1] - q[1], p[2] - q[2] ) ); # Distance between points
PS_ord_y( Q, l ) = (Q.P[ Q.ord_y[ l ] ])[ 2 ];

struct SortByCoord <: Base.Order.Ordering
    P::Vector{Point2F};
    coord::Int
end

import Base.Order.lt
function  lt(o::SortByCoord, a, b)
    p,q = o.P[ a ], o.P[ b ];
    return  (p[o.coord] < q[o.coord] );
end

function  PointSetXY( _P::Points2F )
    PS = PointSetXY( _P, collect( 1:length( _P ) ), collect( 1:length( _P ) ) );

    # Sort the points by x-axis, and store the ordering in PS.ord_x
    sort!( PS.ord_x, order=SortByCoord( PS.P, 1 ) );

    # Sort the points by y-axis, and store the ordering in PS.ord_y
    sort!( PS.ord_y, order=SortByCoord( PS.P, 2 ) );

    return  PS;
end

# Computing subset of points in _PS that fulfill the condition f.
function  PS_filter( P::PointSetXY, f )
    return    PointSetXY( P.P, filter( i -> f( P.P[ i ] ), P.ord_x ),
                             filter( i -> f( P.P[ i ] ), P.ord_y ) );
end

# Use the elevator algorithm to discover any closest pair that is
# better than the current solution.
function  nn_middle( VPS::PointSetXY, mid, sol::Solution )
    low, hi = mid - sol.d, mid + sol.d
    P = PS_filter( VPS, p -> ( low <= p[1] <= hi ) );
    ( len( P ) < 2 )  &&  return  sol;

    b = t = 1;  # Bottom/top of elevator
    while  ( t <= len( P ) )
        # Delete bottom point ∈ elevator if irrelevant
        while  ( b < t ) && ( PS_ord_y( P, b ) < (PS_ord_y( P, t ) - sol.d) )
            b = b + 1;
        end

        # Find nn in points in the elevator to the top point in it.
        # The elevator points are all points with indices in
        # PS.ord_y[ b:t ]
        i_p = P.ord_y[ t ];
        p = P.P[ i_p ];

        for  j ∈ b:(t-1)
            i_q = P.ord_y[ j ]
            q = P.P[ i_q ]
            ℓ = Dist( p, q )
            if  ( ℓ < sol.d )
#                println( "XXXYY" );
                Sol_set( sol, ℓ, i_p, i_q );
            end
        end
        t = t + 1;
    end

    return  sol;
end

function  closest_pair_dc_inner( PS::PointSetXY, sol )
    ( len( PS ) < 2)  &&  return  sol;

    mid = PS.ord_x[ len( PS )  ÷  2 ];
    x_mid = PS.P[ mid ][ 1 ];

    # A subtlety: All the points with x-coordinat x_mid, would be
    # checked by the call nn_middle, so no need to include them in the
    # two recursive calls.

    # Dividing...
    PS_L = PS_filter( PS, p -> ( p[1] < x_mid ) );  # Left points
    PS_R = PS_filter( PS, p -> ( p[1] > x_mid ) );  # Right points

    # Recursing
    closest_pair_dc_inner( PS_L, sol ); # recurse left
    closest_pair_dc_inner( PS_R, sol ); # recurse right

    # Conquering
    nn_middle( PS, x_mid, sol );

    # declaring victory, and withdrawing
    return  sol;
end

function  closest_pair_dc( _PS::Points2F )
    PS = PointSetXY( _PS )

    sol = Solution( Dist(PS.P[1], PS.P[2]), 1, 2 );
    return  closest_pair_dc_inner( PS, sol );
end


############################################################################
# Closest pair in the plane using divide-and-conquer in O( n log n ) time.
############################################################################

function  closest_pair_brute_force( PS::Points2F )
    sol = Solution( Dist(PS[1], PS[2]), 1, 2 );

    for  i in 1:length( PS ) -1
        p = PS[ i ];
        for  j in i+1:length( PS )
            q = PS[ j ];
            ℓ =  Dist( p, q )
            if  ( ℓ < sol.d )
                Sol_set( sol,  ℓ, i, j );
            end
        end
    end

    return  sol;
end

function  main_dc()
    # The D&C algorithm seems faster than the naive algorithm around n ≈ 200
    n = 200;

    println( "n : ", n );
    PS = PointSetXY(  [ (rand(), rand()) for i ∈ 1:n] );

    @time sol = closest_pair( PS );
    @time sol_bf = closest_pair_brute_force( PS );

    println( "cp    : ", sol );
    println( "cp_bf : ", sol_bf );
end

#main();


function  append_to_file( file_path, str )
    try
        open(file_path, "a") do file
            # It's good practice to add a newline character if the appended
            # content should start on a new line.
            write(file,  str)
        end
    catch e
        println("An error occurred while appending to the file: ", e)
    end
end

############################################################################
############################################################################


function  force_compile()
    D=2
    n = 1000
    P = Points( Polygon_random( D,Float64, n ) );
    t_OH_rand = @timed OH_G = OH_closest_pair( P );
    t_dc      = @timed sol = closest_pair_dc( P );
    t_bf      = @timed sol_bf = closest_pair_brute_force( P );
    t_rand    = @timed G = closest_pair( P );
    t_CL_rand = @timed OH_G   = CL_closest_pair( P );
end


function (@main)(ARGS)
    if  length( ARGS ) != 1
        println( "Usage:\n\t"
                 * "closest_pair.jl [n]\n\n" );
        return  -1;
    end

    D = 2;
    n =  parse(Int64, replace( ARGS[ 1  ], "," => "" ) );

    println( "n: ", n );
    force_compile();

    P = Polygon_random( D,Float64, n );

    PB = deepcopy( Points( P ) );
    PC = deepcopy( Points( P ) );
    PD = deepcopy( Points( P ) );
    PE = deepcopy( Points( P ) );
    PF = deepcopy( Points( P ) );

    t_CL_rand = @timed CL_G   = CL_closest_pair( PB );
    t_OH_rand = @timed OH_G   = OH_closest_pair( PC );
    t_dc      = @timed sol    = closest_pair_dc( PD );
    #t_bf      = @timed sol_bf = closest_pair_brute_force( PE );
    t_rand    = @timed G      = closest_pair( PF );

    @printf( "Time D&C           : %10.5f\n", t_dc.time );
    #@printf( "Time Brute force   : %10.5f\n", t_bf.time );
    @printf( "Time Rand          : %10.5f\n", t_rand.time );
    @printf( "Time Rand (OH)     : %10.5f\n", t_OH_rand.time );
    @printf( "Time Rand (CL)     : %10.5f\n", t_CL_rand.time );
    println( "Distance                            : ", sol.d );

    str = @sprintf( "%12d,  %11.6f, %11.6f, %11.6f, %11.6f\n",
        length( P ),
        t_dc.time, t_rand.time,
             t_OH_rand.time, t_CL_rand.time );

    append_to_file( "results/closest_pair.txt", str );

    if  ( sol.d != G.cp_dist )
        println( "\n\n\nBUG!!!!\n\n\n\n" );a
    end
    if  ( sol.d != OH_G.cp_dist )
        println( "\n\n\nBUG in open-hash grid !!!!\n\n\n\n" );
    end
    if  ( sol.d != CL_G.cp_dist )
        println( "\n\n\nBUG in open-hash grid !!!!\n\n\n\n" );
        println( "CL_G.cp_dist: ", CL_G.cp_dist );
    end
    return  0;
end
