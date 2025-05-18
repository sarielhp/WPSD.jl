####################################################
# BBT - BoundingBoxTree
#
# A tree that stores points in a hierarchical axis aligned bounding boxes.
####################################################

module  BBT

push!(LOAD_PATH, pwd()*"/cg/")
push!(LOAD_PATH, pwd() )

using FrechetDist;
using FrechetDist.cg;
using FrechetDist.cg.polygon;
using FrechetDist.cg.point;

#include("./VirtArray.jl")
#include( "VirtArray.jl" );
using VirtArray;
using Cairo, Colors

#using Base: reindex
IntRange = UnitRange{Int};

mutable struct BBTNode{D,T}
    bb::BBox{D,T}
    r::IntRange
    left::Union{Nothing, BBTNode}
    right::Union{Nothing, BBTNode}
    f_leaf::Bool;
    diam::T;
    id::Int
end

mutable struct BBTree{D,T}
    PS::VArray{Point{D,T}};
    root::Union{Nothing, BBTNode}
    id_counter::Int
end

########################################################################

function  node_init( tree::BBTree{D,T}, range::IntRange )  where {D,T}
    bb = BBox{D,T}();
    tree.id_counter += 1;

    #println( "NID: ", tree.id_counter );
    BBox_bound( bb, tree.PS[ range ] );
    return  BBTNode( bb, range, nothing, nothing, false, BBox_diam( bb ),
                     tree.id_counter );
end

function  BBTree_init_inner( p::Polygon{D,T} )::BBTree{D,T} where {D,T}
    varr = VArray( Points( p ) );
    tree = BBTree( varr, nothing, 0 );

    tree.root = node_init( tree, 1:length(p) );

    return  tree;
end

function  BBTree_init_inner( _p::Vector{Point{D,T}} )::BBTree{D,T} where {D,T}
    p = Polygon{D,T}( _p );
    return  BBTree_init_inner( p );
end



function  is_identical( arr::VArray{T}, r )::Bool where {T}
    if  ( length( r ) <= 1 )   return  true;   end

    s = first( r )
    for  i  in  s+1:last(r)
        if  arr[ s ] != arr[ i ]
            return  false;
        end
    end

    return  true;
end


function pnt_varray_partition!( P::VArray{Point{D,T}}, r::IntRange,
                                dim, cutoff )  where {D,T}
    @assert( ( 1 <= first( r ) )  &&  ( last(r) <= length( P ) ) );

    len = top = last( r );
    @inbounds for  i  in  r
        (i > top)  &&  break;
        while  ( P[ i ][ dim ] > cutoff )
            swap!( P, i, top );
            #P[ i ], P[ top ] = P[ top ], P[ i ];
            top = top - 1;
            (i > top)  &&  break;
        end
    end
    return  first(r):top, (top + 1):len;
end


function  node_split( node::BBTNode{D,T}, tree::BBTree{D,T} )  where {D,T}
    if  ( ( ! isnothing( node.left ) )  &&  ( ! isnothing( node.right ) ) )
        #println( "Already split?" );
        return;
    end

    # If a leaf there is nothing to split.
    if  node.f_leaf  return  end


    #println( "node_split" );
    # covers the case range isa 1
    if  is_identical( tree.PS, node.r )
        #println( "Leaf???" );
        node.f_leaf = true;
        return;
    end

    #println( "Doing split!" );
    w, dim = findmax( i->BBox_width( node.bb, i ), 1:D );

    cutoff = BBox_middle( node.bb, dim );
    r_l, r_r = pnt_varray_partition!( tree.PS, node.r, dim, cutoff )
    @assert( ( length( r_l ) > 0 )  &&  ( length( r_r ) > 0 ) );

    node.left = node_init( tree, r_l );
    node.right = node_init( tree, r_r );
    #println( "L NID: ", node.left.id )
    #println( "R NID: ", node.right.id )
    
    return  node;
end


function  node_expand( v::BBTNode{D,T}, tree::BBTree{D,T} ) where {D,T}
    v.f_leaf  &&  return;

    node_split( v, tree );

    #    v.f_leaf  &&  return;
    #
#    ( v.left != nothing )   &&  node_expand( v.left, tree );
#    ( v.right != nothing )  &&  node_expand( v.right, tree );
end

function BBTree_init( PS::Polygon{D,T} )::BBTree{D,T} where {D,T}
    tree = BBTree_init_inner( PS );

    #println( typeof( tree ) );
    #println( "tree created" );

    node_expand( tree.root, tree );

    return  tree;
end

function  fully_expand( node::BBTNode{D,T}, tree::BBTree{D,T} ) where {D,T}
    if   node.f_leaf  return  end;

    if ( isnothing( node.left )  ||  isnothing( node.right ) )
        #println( "   mode expand? " );
        node_expand( node, tree );
    end
    if  ( ! isnothing( node.left ) )
        fully_expand( node.left, tree );
    end
    if  ( ! isnothing( node.right ) )
        fully_expand( node.right, tree );
    end
end

function  BBTree_fully_expand( tree::BBTree{D,T} ) where {D,T}
    fully_expand( tree.root, tree );
end

function   bbox_draw( context, bb, color )
    set_source(context, color)

    bl = BBox_bottom_left( bb );
    w = BBox_width( bb );
    h = BBox_height( bb );

    rectangle(context, bl[1], bl[2], w, h );
    fill( context );
end


function   node_draw( context, node, level, range )
    if  node.f_leaf  return  end;

    if  level ∈ range
        # 0.1 alpha means 10% opaque, 90% transparent
        yellow_transparent = coloralpha(colorant"yellow", 0.1)

        set_source(context, yellow_transparent)

        bb = BBox_expand( node.bb, 1.01 );

        bl = BBox_bottom_left( bb );
        w = BBox_width( bb );
        h = BBox_height( bb );

        # Draw the rectangle
        rectangle(context, bl[1], bl[2], w, h );
        fill_preserve(context)

        # Draw the blue frame
        set_source(context, colorant"blue")
        set_line_width(context, 1.0) # Set the line width for the frame
        stroke(context)
    end
    if  ( ! isnothing( node.left ) )
        node_draw( context, node.left, level+1, range );
    end
    if  ( ! isnothing( node.right ) )
        node_draw( context, node.right, level+1, range );
    end
end

function  depth( node::BBTNode{D,T} )::Int where  {D,T}
    if  ( isnothing( node ) )  return  0  end;
    if  node.f_leaf  return 1 end;

    return 1 + max(  depth( node.left ), depth( node.right ) )
end

function BBTree_draw( tree::BBTree{D,T}, filename::String ) where{D,T}
    bb = BBox_expand( tree.root.bb, 1.3 );
    #BBox_print( bb );

    surface = CairoPDFSurface(filename, BBox_max( bb, 1 ),
        BBox_max( bb, 2 ) )
    context = CairoContext(surface)

    bbox_draw( context, bb, colorant"lightblue" );

    d = depth( tree.root );
    node_draw( context, tree.root, 0, 0:d );

    # Set the fill color to yellow with 90% transparency

    show_page(context)

    # Now draw levels of the tree...
    for  i ∈ (d-1):-1:0
        bbox_draw( context, bb, colorant"lightblue" );
        node_draw( context, tree.root, 0, i:i+1 );
        show_page( context )
    end

    finish(surface)
    println( "\n" * "Created file: $filename")
end


function  BBTree_refine_node( tree, node )
    node_split( node, tree );
end

# Export list...

export  BBTNode,  BBTree;
export  BBTree_draw, BBTree_init, BBTree_fully_expand
export  BBTree_refine_node


end  # Module
####################################################################3
