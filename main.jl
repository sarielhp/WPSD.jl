####################################################
# main: Test all the new modules...
#
####################################################


push!(LOAD_PATH, pwd()*"/cg/")
push!(LOAD_PATH, pwd() )

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
using VirtArray;
using Printf;

include( "BBT.jl" );
include( "WSPD.jl" );

using .BBT
using .WSPD
using Cairo, Colors

function  BBT_test( P, N )
    tree = BBTree_init( P );
    #println( "Tree initialized..." );
    BBTree_fully_expand( tree );
    BBTree_draw( tree, "test.pdf" );
end

function  exact_diameter( P::Polygon{D,T}  ) where {D,T}
    curr::Float64 = 0;
    
    for  i ∈ 1:length(P)-1
        for j ∈ i+1:length(P)
            curr = max( curr, Dist( P[i], P[j] ) );
        end
    end
    return  curr;
end


function  approx_diameter( P::Polygon{D,T}, ε::Float64  ) where {D,T}
    if  length( P ) <= 1
        return  0.0
    end
    c = 1.0 + ε;
    W = WSPD_init( P,  ε / 2  );
    curr = Dist( P[ 1 ], P[ 2 ] );
    while  ( WSPD_top_diam_ub( W ) > c * curr )
        top = WSPD_get_top( W );
        p,q = WSPD_get_reps( W, top );
        curr = max( curr, Dist( p, q ) );
        WSPD_top_refine( W );
    end

    return  curr
end

function  test_diameter( N::Int64, f_silent = false )
    P = Polygon_random( 2, Float64, N );
    mult!( P, 900.0 );
    shift!( P, Point( 200.0, 200.0 ) );

#    println( "----------------------------------------------" );
#    println( "number of points ", N );
#    println( "Approx: " );
    t_approx = @timed approx_diam = approx_diameter( P, 0.001 );

#    println( "Exact: " );
    t_exact = @timed exact_diam = exact_diameter( P );

    @printf( "N: %8d  diam ≈ %10g  Exact: %10g    T≈ %10.4f  ET= %10.4f\n",
             N, approx_diam, exact_diam, t_approx.time, t_exact.time );
#    println( "N: ", N, " diam ≈ ", approx_diam, "  exact: ", exact_diam );
#    println( "exact  diam: ", exact_diam );
end

function (@main)(ARGS)
    for  i ∈ 1:20
        N = 2^i 
        test_diameter( N );
    end
    exit( -1 );
    
    N = 400000;
    P = Polygon_random( 2, Float64, N );

    mult!( P, 900.0 );
    shift!( P, Point( 200.0, 200.0 ) );

    println( "N: ", N );
    #W = WSPD_init( P,  0.1 );
    #WSPD_expand( W );
    println( "Approx: " );
    t = @time approx_diam = approx_diameter( P, 0.1 );
    println( "TTT: ", t );
    println( "approx diam: ", approx_diam );

    println( "Exact: " );
    @time exact_diam = exact_diameter( P );
    println( "exact  diam: ", exact_diam );


    #println( "Pairs computed: ", length( W.finals ) );
    #println( "Pairs computed: ", length( W.finals )/ (N*N) );
    return  0;
end


####################################################################3

#export  BBTree_build;

#end
