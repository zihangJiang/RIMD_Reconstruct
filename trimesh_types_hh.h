#ifndef TRIMESH_TYPES_HH_H
#define TRIMESH_TYPES_HH_H

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
//== TYPEDEFS =================================================================

/** Default traits for the PolyMesh
*/
struct TriTraits : public OpenMesh::DefaultTraits
{
  /// Use double precision points
  typedef OpenMesh::Vec3d Point;
  /// Use double precision Normals
  typedef OpenMesh::Vec3d Normal;
  /// Use double precision TexCood2D
  typedef OpenMesh::Vec2d TexCoord2D;

  /// Use RGBA Color
  typedef OpenMesh::Vec4f Color;
};

/// Simple Name for Mesh
typedef OpenMesh::TriMesh_ArrayKernelT<TriTraits>  TriMesh;
#endif // TRIMESH_TYPES_HH_H
