#include <Deformable/GaussLegendre.h>
#include <Deformable/Hexahedron.h>

USE_PRJ_NAMESPACE

int main()
{
  std::function<scalarD(Vec3d)> f=[&](Vec3d) {
    return 1;
  };
  for(sizeType it=0; it<10; it++) {
    Vec3d a=Vec3d::Random();
    Vec3d b=Vec3d::Random();
    Vec3d c=Vec3d::Random();
    scalarD val1=GaussLegendreIntegral<scalarD,scalarD>::integrateTri(a,b,c,f,1);
    scalarD val2=GaussLegendreIntegral<scalarD,scalarD>::integrateTri(a,b,c,f,2);
    scalarD val3=GaussLegendreIntegral<scalarD,scalarD>::integrateTri(a,b,c,f,3);
    INFOV("val(deg)=%f(1),%f(2),%f(3) area(tri)=%f",val1,val2,val3,TriangleTpl<scalarD>(a,b,c).area())
  }
  for(sizeType it=0; it<10; it++) {
    Vec3d a=Vec3d::Random();
    Vec3d b=Vec3d::Random();
    Vec3d c=Vec3d::Random();
    Vec3d d=Vec3d::Random();
    scalarD val1=GaussLegendreIntegral<scalarD,scalarD>::integrateTet(a,b,c,d,f,1);
    scalarD val2=GaussLegendreIntegral<scalarD,scalarD>::integrateTet(a,b,c,d,f,2);
    scalarD val3=GaussLegendreIntegral<scalarD,scalarD>::integrateTet(a,b,c,d,f,3);
    INFOV("val(deg)=%f(1),%f(2),%f(3) volume(tet)=%f",val1,val2,val3,TetrahedronTpl<scalarD>(a,b,c,d).volume())
  }

  //hexahedron
  HexahedronTpl<scalarD> hex(Vec3d(0,0,0),Vec3d(1,0,0),Vec3d(0,1,0),Vec3d(1,1,0),
                             Vec3d(0,0,1),Vec3d(1,0,1),Vec3d(0,1,1),Vec3d(1,1,1));
  scalarD val=GaussLegendreIntegral<scalarD,scalarD>::integrateHex
              (Vec3d(0,0,0),Vec3d(1,0,0),Vec3d(0,1,0),Vec3d(1,1,0),
               Vec3d(0,0,1),Vec3d(1,0,1),Vec3d(0,1,1),Vec3d(1,1,1),
               f,1);
  INFOV("volumeRef(hex)=%f, volume(hex)=%f",val,hex.volume())
  HexahedronTpl<scalarD> hex2(Vec3d(0,0,0)+Vec3d::Random()*0.1f,
                              Vec3d(1,0,0)+Vec3d::Random()*0.1f,
                              Vec3d(0,1,0)+Vec3d::Random()*0.1f,
                              Vec3d(1,1,0)+Vec3d::Random()*0.1f,
                              Vec3d(0,0,1)+Vec3d::Random()*0.1f,
                              Vec3d(1,0,1)+Vec3d::Random()*0.1f,
                              Vec3d(0,1,1)+Vec3d::Random()*0.1f,
                              Vec3d(1,1,1)+Vec3d::Random()*0.1f);
  hex2.debugDistanceVTK(1,30,"hexDist.vtk");
  hex2.writeVTK("hex.vtk");
  return 0;
}
