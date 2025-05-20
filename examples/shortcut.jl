#! /usr/bin/env  julia

# shortcut.jl:
#
# Approximate the optimal shortcut for a give data file.
#

push!(LOAD_PATH, pwd()*"/src/")
push!(LOAD_PATH, pwd()*"/src/cg/")

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;
using VirtArray;
using Printf;
using DataFrames
using PrettyTables

include( "graphics.jl" )

using BBT
using WSPD
using Cairo, Colors


function  rt_str( rt )
    return  @sprintf( "%.6f", rt );
end

function  rt_str_2( rt )
    return  @sprintf( "%.2f", rt );
end

function  approx_shortcut( P::Polygon{D,T}, ε::Float64  ) where {D,T}
    @assert( length( P ) > 1 );
    
    if  length( P ) <= 1
        return  0.0
    end

    plens = Polygon_prefix_lengths( P );

    c = 1.0 + ε;
    W = WSPD_init( P,  ε / 2  );
    WSPD_expand( W );

    max_quality = -1.0;
    p_i = q_i = 1;
    for  pair ∈ W.finals
        l_i,r_i = WSPD_reps_orig_indexes( W, pair );
        l = abs( plens[ l_i ] - plens[ r_i ] )
        d = Dist( P[ l_i ], P[ r_i ] );

        if  d == 0.0   continue  end

        quality = l / d
        if quality > max_quality
            p_i = min( l_i, r_i);
            q_i = max( l_i, r_i );
            max_quality = quality
        end
    end
    
    return  p_i, q_i, max_quality;
end

function (@main)(ARGS)
    if  length( ARGS ) != 2
        println( "Usage:\n\t"
                 * "shortcut.jl [eps] [file.txt]\n\n" );
        return  -1;
    end
#    println( ARGS[ 1 ] );
    eps =  parse(Float64, ARGS[ 1  ] );
    P = read_file( ARGS[ 2 ] );

    P = Polygon_sample_uniformly( P, 1000 );
    println( "N = ", length( P ) );
    approx_shortcut( P, 1.0 );
    list = VecPolygon2F();
    for  i ∈ 1:40
        push!( list, P );
        t = @timed  v_i,v_j,quality = approx_shortcut( P, eps );
        println( i, " : ", quality );
        if  ( quality < 1.250 )
            break;
        end
        println( v_i, "..", v_j, " :   Time: ", t.time );
        P = Polygon_shortcut( P, v_i, v_j );
    end
    
    #println( i, "   ", j );
    #println( "Quality: ", quality );

    #shortcut = Polygon2F()
    #push!( shortcut, P[ i ], P[j] );

    #push!( list, P, shortcut );

    bb = BBox2F();
    BBox_bound( bb, P );


    p = BBox_bottom_left( bb );
    for  poly  in list
        Polygon_translate!( poly, p );
    end

    output_polygons_to_file( list, "curves.pdf", true );

    
#    @printf( "%.6f-approx shortcut: %15.6f   Time: %11.6f    N: %7d\n",
#             eps, diam, t.time, length( P ) );
#    println( eps"-approx diameter: ", diam,  " time: " );
    return  0;
end


####################################################################3

#export  BBTree_build;

#end
