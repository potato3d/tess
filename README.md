# tess
Tessellation routines for parametric 3D geometries

# Description

The main files tessellator.* generate optimized 3D triangle meshes for the following geometries:
* box
* circular_torus
* cone
* cone_offset
* cylinder
* cylinder_offset
* dish
* pyramid
* rectangular_torus
* sphere

Checkout geometries.h for details on how each geometry can be parameterized.

Polygon tessellation is implemented in the polygon_tessellator.* files. It uses the GLUT tessellation library.

The mesh_builder.* files convert a triangle mesh represented as GL_TRIANGLE_FAN or GL_TRIANGLE_STRIP into GL_TRIANGLES representation.

The mesh_optimizer.* files implement algorithms to optimize an existing triangle mesh. It can remove unused vertices in O(n) and merge vertices with exactly the same coordinates in O(n) (using a hash table).