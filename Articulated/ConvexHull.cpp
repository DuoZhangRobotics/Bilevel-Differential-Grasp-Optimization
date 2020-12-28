#ifdef CGAL_SUPPORT
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/print_wavefront.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <fstream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef K::Point_3 Point_3;
#include "ConvexHull.h"

PRJ_BEGIN

void readOFF(std::istream& is,ObjMeshD& mesh)
{
  Polyhedron poly;
  is >> poly;
  {
    std::ofstream os("tmp.obj");
    print_polyhedron_wavefront(os,poly);
  }
  std::ifstream isObj("tmp.obj");
  mesh.read(isObj,false,false);
}
//convexification
ObjMeshF makeConvex(const ObjMeshF& in)
{
  ObjMeshD inD;
  ObjMeshF out;
  in.cast<scalarD>(inD);
  makeConvex(inD).cast<scalarF>(out);
  return out;
}
ObjMeshD makeConvex(const ObjMeshD& in)
{
  std::vector<Point_3> points;
  for(sizeType i=0; i<(sizeType)in.getV().size(); i++) {
    Point_3 pt(in.getV(i)[0],in.getV(i)[1],in.getV(i)[2]);
    points.push_back(pt);
  }
  Polyhedron_3 poly;
  CGAL::convex_hull_3(points.begin(),points.end(),poly);
  {
    std::ofstream os("tmp.off");
    os << poly;
  }
  ObjMeshD out;
  std::ifstream is("tmp.off");
  readOFF(is,out);

  out.smooth();
  out.makeUniform();
  if(out.getVolume()<0)
    out.insideOut();
  return out;
}

PRJ_END
#else
#include <CommonFile/ObjMesh.h>
PRJ_BEGIN

void readOFF(std::istream& is,ObjMeshD& mesh)
{
  FUNCTION_NOT_IMPLEMENTED
}
//convexification
ObjMeshF makeConvex(const ObjMeshF& in)
{
  FUNCTION_NOT_IMPLEMENTED
}
ObjMeshD makeConvex(const ObjMeshD& in)
{
  FUNCTION_NOT_IMPLEMENTED
}

PRJ_END
#endif
