#! /usr/bin/julia

push!(LOAD_PATH, pwd() )
#using Printf;
#using Graphs;

include( "cg.jl" );
using .cg;

include( "graphics.jl" );
include( "BBT.jl" );
#using .BBT

function  main( N )
    PS = Polygon_random( 2, Float64, N );
    tree = BBTree_build( PS );
end

main( 10 );

