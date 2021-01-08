#include <Articulated/ArticulatedBody.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

#define NORMAL Vec3(0,0,1)
#define TANGENT1 Vec3(1,0,0)
#define TANGENT2 Vec3(0,1,0)
#define SCALE 1.0
#define SIZE 10
void drawFloor(const std::string& path,scalar LEVEL)
{
  create(path);
  for(sizeType i=0; i<2; i++) {
    std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
    VTKWriter<scalar> os("floor",path+"/floor"+std::to_string(i)+".vtk",true);
    sizeType idx=0;
    for(scalar x=-SIZE; x<=SIZE-SCALE; x+=SCALE,idx++) {
      sizeType idy=0;
      for(scalar y=-SIZE; y<=SIZE-SCALE; y+=SCALE,idy++)
        if((idx+idy)%2==i) {
          vss.push_back(TANGENT1*x+TANGENT2*y+NORMAL*(LEVEL-SCALE));
          vss.push_back(TANGENT1*x+TANGENT2*y+NORMAL*(LEVEL-SCALE)+Vec3::Constant(SCALE));
        }
    }
    os.appendVoxels(vss.begin(),vss.end(),true);
  }
}
int main()
{
  drawFloor("Robosimian",-0.8f);
  drawFloor("Spider",-0.435);
  drawFloor("Bird",-0.81);
  return 0;
}
