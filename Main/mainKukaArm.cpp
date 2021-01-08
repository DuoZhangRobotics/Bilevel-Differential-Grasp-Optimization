#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Main/GenericArticulatedObjective.h>
#include <Utils/DebugGradient.h>

USE_PRJ_NAMESPACE

typedef double T;
int main()
{
  //read the robot
  typedef PBDArticulatedGradientInfo<T>::Vec Vec;
  typedef PBDArticulatedGradientInfo<T>::Vec3T Vec3T;
  typedef PBDArticulatedGradientInfo<T>::Mat3XT Mat3XT;
  typedef PBDArticulatedGradientInfo<T>::Mat3X4T Mat3X4T;
  ArticulatedBody body=ArticulatedLoader::readURDF("data/kuka_lwr/kuka.urdf",false,true);

  //compute robot pose and draw it
  Vec theta=Vec::Random(body.nrDOF());  //robot parameter (theta)
  PBDArticulatedGradientInfo<T> info(body,theta);
  body.writeVTK(info._TM,"pose.vtk",Joint::MESH);

  //get joint limit:
  std::cout << "LowerLimit: " << body.lowerLimit().transpose() << std::endl;
  std::cout << "UpperLimit: " << body.upperLimit().transpose() << std::endl;

  //get the mesh of a certain body link
  ObjMesh m=body.joint(0).getMesh(Joint::MESH);
  m.writeVTK("joint0.vtk",true);

  //compute the global position a certain vertex on body link
  sizeType JID=2;   //which joint are you considering
  m=body.joint(JID).getMesh(Joint::MESH);
  Vec3T vLocal=m.getV(0).template cast<T>();
  Mat3X4T trans=TRANSI(info._TM,JID);   //the 3x4 local->global joint transformation matrix
  Vec3T vGlobal=ROT(trans)*vLocal+CTR(trans);
  std::cout << "vertex in local coordinates: " << vLocal.transpose() << std::endl;
  std::cout << "vertex in world coordinates: " << vGlobal.transpose() << std::endl;

  {
    //compute the derivatives of global position with respect to theta
    Mat3XT J=Mat3XT::Zero(3,theta.size());
    info.JRCSparse(body,JID,[&](sizeType col,const Vec3T& w) {
      J.col(col)+=w.cross(ROTI(info._TM,JID)*vLocal);
    },[&](sizeType col,const Vec3T& w) {
      J.col(col)+=w;
    });
    //compare the derivative with finite difference
    DEFINE_NUMERIC_DELTA_T(T)
    Vec Dtheta=Vec::Random(body.nrDOF());
    PBDArticulatedGradientInfo<T> info2(body,theta+Dtheta*DELTA);
    Mat3X4T trans2=TRANSI(info2._TM,JID);
    Vec3T vGlobal2=ROT(trans2)*vLocal+CTR(trans2);
    DEBUG_GRADIENT("DvDtheta",std::sqrt((J*Dtheta).squaredNorm()),std::sqrt((J*Dtheta-(vGlobal2-vGlobal)/DELTA).squaredNorm()))
  }

  {
    //now we illustrate the hessian using a template function GenericArticulatedEnergy:
    ExampleArticulatedObjective<T>(body,100).debug(theta);
  }
  return 0;
}
