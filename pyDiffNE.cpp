#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <TrajOpt/Environment/EnvWrenchConstructor.h>
#include <TrajOpt/Environment/GranularWrenchConstructor.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedUtils.h>
#include <Articulated/MDPSimulator.h>
#include <Articulated/RotationUtil.h>
#include <Articulated/PDTarget.h>
#include <CommonFile/geom/StaticGeom.h>
#include <sstream>

USE_PRJ_NAMESPACE
namespace py=pybind11;
PYBIND11_DECLARE_HOLDER_TYPE(T,std::shared_ptr<T>);
#define PYTHON_FUNC(C,NAME) .def(STRINGIFY_OMP(NAME),&C::NAME)
#define PYTHON_SFUNC(C,NAME) .def_static(STRINGIFY_OMP(NAME),&C::NAME)
#define PYTHON_OFUNC(C,NAME,SIG) .def(STRINGIFY_OMP(NAME),SIG &C::NAME)
#define PYTHON_OSFUNC(C,NAME,SIG) .def_static(STRINGIFY_OMP(NAME),SIG &C::NAME)
#define PYTHON_MEMBER_READ(C,NAME) .def_readonly(STRINGIFY_OMP(NAME),&C::NAME)
#define PYTHON_MEMBER_READWRITE(C,NAME) .def_readwrite(STRINGIFY_OMP(NAME),&C::NAME,&C::NAME)
#define PYTHON_FIELD_READWRITE(C,NAME) .def_readwrite(STRINGIFY_OMP(NAME),&C::NAME)

#define EIGEN_OP_VEC(T) \
.def(py::init<>())  \
.def(py::init([](const T& other) {return T(other);}))  \
.def(py::init([](py::list lst) {T v;v.resize(lst.size());sizeType id=0;for(auto item:lst)v[id++]=item.cast<T::Scalar>();return v;}))  \
.def("__deepcopy__",[](const T& v,py::dict){return T(v);})   \
.def("__add__",[](const T& v,const T& w){return T(v+w);},py::is_operator()) \
.def("__sub__",[](const T& v,const T& w){return T(v-w);},py::is_operator()) \
.def("__mul__",[](const T& v,const T::Scalar& c){return T(v*c);},py::is_operator()) \
.def("__div__",[](const T& v,const T::Scalar& c){return T(v/c);},py::is_operator()) \
.def("fromList",[](T& v,py::list lst){v.resize(lst.size());sizeType id=0;for(auto item:lst)v[id++]=item.cast<T::Scalar>();}) \
.def("toList",[](const T& v){py::list li;for(sizeType i=0;i<v.size();i++)li.append(v[i]);return li;})  \
.def("castFrom",[](T& v,const T& v2){v=v2.template cast<T::Scalar>();}) \
.def("castToDouble",[](const T& v){return v.template cast<double>().eval();}) \
.def("setZero",[](T& v){v.setZero();})  \
.def("setOnes",[](T& v){v.setOnes();})  \
.def("setConstant",[](T& v,T::Scalar val){v.setConstant(val);})  \
.def("resize",[](T& v,sizeType sz){v.resize(sz);})   \
.def("dot",[](const T& v,const T& w){return v.dot(w);}) \
.def("size",[](const T& v){return v.size();})   \
.def("__len__",[](const T& v){return v.size();})   \
.def("__getitem__",[](const T& v,sizeType i){return v[i];}) \
.def("__setitem__",[](T& v,sizeType i,T::Scalar val){v[i]=val;})    \
.def("__repr__",[](const T& v){std::ostringstream os;os<<v.unaryExpr([&](const typename T::Scalar& in) {return (scalar)std::to_double(in);});return os.str();});

#define EIGEN_OP_MAT(T) \
.def(py::init<>())  \
.def(py::init([](const T& other) {return T(other);}))  \
.def(py::init([](py::list lst) {T v;v.resize(lst.size(),lst.begin()->template cast<py::list>().size());sizeType id=0;for(auto item:lst){sizeType id2=0;for(auto item2:item)v(id,id2++)=item2.template cast<T::Scalar>();id++;}return v;}))  \
.def("__deepcopy__",[](const T& v,py::dict){return T(v);})   \
.def("__add__",[](const T& v,const T& w){return T(v+w);},py::is_operator()) \
.def("__sub__",[](const T& v,const T& w){return T(v-w);},py::is_operator()) \
.def("__mul__",[](const T& v,const T::Scalar& c){return T(v*c);},py::is_operator()) \
.def("__div__",[](const T& v,const T::Scalar& c){return T(v/c);},py::is_operator()) \
.def("fromList",[](T& v,py::list lst){v.resize(lst.size(),lst.begin()->template cast<py::list>().size());sizeType id=0;for(auto item:lst){sizeType id2=0;for(auto item2:item)v(id,id2++)=item2.template cast<T::Scalar>();id++;}}) \
.def("toList",[](const T& v){py::list li;for(sizeType i=0;i<v.rows();i++){py::list lic;for(sizeType j=0;j<v.cols();j++)lic.append(v(i,j));li.append(lic);}return li;})  \
.def("castFrom",[](T& v,const T& v2){v=v2.template cast<T::Scalar>();}) \
.def("castToDouble",[](const T& v){return v.template cast<double>().eval();}) \
.def("setZero",[](T& v){v.setZero();})  \
.def("setOnes",[](T& v){v.setOnes();})  \
.def("setConstant",[](T& v,T::Scalar val){v.setConstant(val);})  \
.def("resize",[](T& v,sizeType r,sizeType c){v.resize(r,c);})   \
.def("rows",[](const T& v){return v.rows();})   \
.def("cols",[](const T& v){return v.cols();})   \
.def("__len__",[](const T& v){return v.size();})   \
.def("__getitem__",[](const T& v,sizeType i,sizeType j){return v(i,j);}) \
.def("__getitem__",[](const T& v,py::tuple t){return v(t[0].template cast<sizeType>(),t[1].template cast<sizeType>());}) \
.def("__setitem__",[](T& v,sizeType i,sizeType j,T::Scalar val){v(i,j)=val;})    \
.def("__setitem__",[](T& v,py::tuple t,T::Scalar val){v(t[0].template cast<sizeType>(),t[1].template cast<sizeType>())=val;})    \
.def("__repr__",[](const T& v){std::ostringstream os;os<<v.unaryExpr([&](const typename T::Scalar& in) {return (scalar)std::to_double(in);});return os.str();});

#define ALL_MAT_VEC(T,SUFFIX)  \
py::class_<MDPSimulator<T>::Vec3T>(m,("Vec3"+std::string(SUFFIX)).c_str())   \
.def(py::init<T,T,T>()) \
EIGEN_OP_VEC(MDPSimulator<T>::Vec3T)    \
py::class_<MDPSimulator<T>::Vec4T>(m,("Vec4"+std::string(SUFFIX)).c_str())   \
.def(py::init<T,T,T,T>())   \
EIGEN_OP_VEC(MDPSimulator<T>::Vec4T)    \
py::class_<MDPSimulator<T>::Vec6T>(m,("Vec6"+std::string(SUFFIX)).c_str())   \
EIGEN_OP_VEC(MDPSimulator<T>::Vec6T)    \
py::class_<MDPSimulator<T>::Vec>(m,("Vec"+std::string(SUFFIX)).c_str())  \
EIGEN_OP_VEC(MDPSimulator<T>::Vec)  \
py::class_<MDPSimulator<T>::Mat3T>(m,("Mat3"+std::string(SUFFIX)).c_str())   \
EIGEN_OP_MAT(MDPSimulator<T>::Mat3T)  \
py::class_<MDPSimulator<T>::Mat3XT>(m,("Mat3X"+std::string(SUFFIX)).c_str())   \
EIGEN_OP_MAT(MDPSimulator<T>::Mat3XT)  \
py::class_<MDPSimulator<T>::Mat3X4T>(m,("Mat3X4"+std::string(SUFFIX)).c_str())   \
EIGEN_OP_MAT(MDPSimulator<T>::Mat3X4T)  \
py::class_<MDPSimulator<T>::Mat6XT>(m,("Mat6"+std::string(SUFFIX)).c_str()) \
EIGEN_OP_MAT(MDPSimulator<T>::Mat6XT)   \
py::class_<MDPSimulator<T>::MatT>(m,("Mat"+std::string(SUFFIX)).c_str()) \
EIGEN_OP_MAT(MDPSimulator<T>::MatT)

namespace pybind11 {
namespace detail {
template <>
struct type_caster<__float128> : public type_caster_base<__float128>
{
  using base=type_caster_base<__float128>;
public:
  bool load(handle src,bool convert) {
    if(base::load(src,convert))
      return true;
    else if(py::isinstance<py::int_>(src)) {
      value=new __float128(py::cast<int>(src));
      return true;
    } else if(py::isinstance<py::float_>(src)) {
      value=new __float128(py::cast<double>(src));
      return true;
    }
    return false;
  }
  static handle cast(const __float128& src,return_value_policy policy,handle parent) {
    return base::cast(std::to_double(src),policy,parent);
  }
};
template <>
struct type_caster<mpfr::mpreal> : public type_caster_base<mpfr::mpreal>
{
  using base=type_caster_base<mpfr::mpreal>;
public:
  bool load(handle src,bool convert) {
    if(base::load(src,convert))
      return true;
    else if(py::isinstance<py::int_>(src)) {
      value=new mpfr::mpreal(py::cast<int>(src));
      return true;
    } else if(py::isinstance<py::float_>(src)) {
      value=new mpfr::mpreal(py::cast<double>(src));
      return true;
    }
    return false;
  }
  static handle cast(const mpfr::mpreal& src,return_value_policy policy,handle parent) {
    return base::cast(std::to_double(src),policy,parent);
  }
};
}
}
PYBIND11_MODULE(pyDiffNE,m)
{
  //-----------------------------------------------------------------basic types
  ALL_MAT_VEC(sizeType,"i")
  ALL_MAT_VEC(float,"f")
  ALL_MAT_VEC(double,"d")
  ALL_MAT_VEC(__float128,"q")
  ALL_MAT_VEC(mpfr::mpreal,"m")
  m.def("eulerX1Y3Z2",[](const Vec3d& w) {
    return eulerX1Y3Z2<scalarD,Vec3d>(w,NULL,NULL);
  });
  m.def("setNumThreads",[](sizeType t) {
    ASSERT_MSG(t>=1 && t<100,"Number of threads must be within: [1,100)!")
    OmpSettings::getOmpSettingsNonConst().setNrThreads(t);
    INFOV("DiffNE using %ld thread!",t)
  });
  //-----------------------------------------------------------------basic numerics
#define MULTI_PREC_SCALAR(T,SUFFIX)   \
py::class_<T>(m,std::string(SUFFIX).c_str()) \
.def(py::init<T>()) \
.def(py::init<double>())    \
.def(py::init<float>()) \
.def(py::init<>())  \
.def("__add__",[](const T& v,const T& w) {return v+w;},py::is_operator())   \
.def("__radd__",[](const T& v,const T& w) {return w+v;},py::is_operator())  \
.def("__sub__",[](const T& v,const T& w) {return v-w;},py::is_operator())   \
.def("__rsub__",[](const T& v,const T& w) {return w-v;},py::is_operator())  \
.def("__mul__",[](const T& v,const T& c) {return v*c;},py::is_operator())   \
.def("__rmul__",[](const T& v,const T& c) {return c*v;},py::is_operator())  \
.def("__div__",[](const T& v,const T& c) {return v/c;},py::is_operator())   \
.def("__rdiv__",[](const T& v,const T& c) {return c/v;},py::is_operator())  \
.def("__add__",[](const T& v,double w) {return v+w;},py::is_operator())     \
.def("__radd__",[](const T& v,double w) {return w+v;},py::is_operator())    \
.def("__sub__",[](const T& v,double w) {return v-w;},py::is_operator())     \
.def("__rsub__",[](const T& v,double w) {return w-v;},py::is_operator())    \
.def("__mul__",[](const T& v,double c) {return v*c;},py::is_operator())     \
.def("__rmul__",[](const T& v,double c) {return c*v;},py::is_operator())    \
.def("__div__",[](const T& v,double c) {return v/c;},py::is_operator())     \
.def("__rdiv__",[](const T& v,double c) {return c/v;},py::is_operator())    \
.def("__repr__",[](const T& v) {std::ostringstream os;os<<std::to_string(v);return os.str();})  \
.def("__float__",[](const T& v) {return std::to_double(v);})    \
.def_static("set_prec",[](unsigned prec) {mpfr_set_default_prec(prec);});
  MULTI_PREC_SCALAR(__float128,"float128")
  MULTI_PREC_SCALAR(mpfr::mpreal,"mpreal")
  //-----------------------------------------------------------------VTKWriter
  py::class_<VTKWriter<scalar>>(m,"VTKWriter").def(py::init<std::string,std::string,bool>());
  //-----------------------------------------------------------------ObjMesh
  py::class_<ObjMesh>(m,"ObjMesh")
  .def(py::init<>())
  PYTHON_FUNC(ObjMesh,read)
  PYTHON_OFUNC(ObjMesh,write,(bool(ObjMesh::*)(const std::string&)const))
  PYTHON_OFUNC(ObjMesh,writeVTK,(bool(ObjMesh::*)(const std::string&,bool,bool,bool,const Vec4*,const std::vector<scalar>*)const))
  PYTHON_FUNC(ObjMesh,getTC)
  PYTHON_OFUNC(ObjMesh,getV,(const Vec3&(ObjMesh::*)(int)const))
  PYTHON_OFUNC(ObjMesh,getI,(const Vec3i&(ObjMesh::*)(int)const))
  PYTHON_OFUNC(ObjMesh,getIG,(const int&(ObjMesh::*)(int)const))
  PYTHON_OFUNC(ObjMesh,getN,(const Vec3&(ObjMesh::*)(int)const))
  PYTHON_OFUNC(ObjMesh,getTN,(const Vec3&(ObjMesh::*)(int)const))
  PYTHON_OFUNC(ObjMesh,getArea,(scalar(ObjMesh::*)(int)const))
  PYTHON_OFUNC(ObjMesh,getV,(const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >&(ObjMesh::*)()const))
  PYTHON_OFUNC(ObjMesh,getI,(const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >&(ObjMesh::*)()const))
  PYTHON_OFUNC(ObjMesh,getIG,(const std::vector<int>&(ObjMesh::*)()const))
  PYTHON_OFUNC(ObjMesh,getN,(const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >&(ObjMesh::*)()const))
  PYTHON_OFUNC(ObjMesh,getFN,(const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >&(ObjMesh::*)()const))
  PYTHON_OFUNC(ObjMesh,getTN,(const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >&(ObjMesh::*)()const))
  PYTHON_OFUNC(ObjMesh,getT,(const Mat3&(ObjMesh::*)()const))
  PYTHON_OFUNC(ObjMesh,getT,(Mat3&(ObjMesh::*)()))
  PYTHON_OFUNC(ObjMesh,getPos,(const Vec3&(ObjMesh::*)()const))
  PYTHON_OFUNC(ObjMesh,getPos,(Vec3&(ObjMesh::*)()))
  PYTHON_OFUNC(ObjMesh,getScale,(const scalar&(ObjMesh::*)()const))
  PYTHON_OFUNC(ObjMesh,getScale,(scalar&(ObjMesh::*)()))
  PYTHON_OFUNC(ObjMesh,applyTrans,(void(ObjMesh::*)()))
  PYTHON_OFUNC(ObjMesh,applyTrans,(void(ObjMesh::*)(const Vec3&)))
  PYTHON_FUNC(ObjMesh,getBB)
  PYTHON_FUNC(ObjMesh,getId)
  PYTHON_FUNC(ObjMesh,getVolume)
  PYTHON_FUNC(ObjMesh,getCentroid)
  PYTHON_FUNC(ObjMesh,getVolumeCentroid)
  PYTHON_OFUNC(ObjMesh,centroidOffset,(Vec3&(ObjMesh::*)()))
  PYTHON_OFUNC(ObjMesh,centroidOffset,(const Vec3&(ObjMesh::*)()const))
  PYTHON_FUNC(ObjMesh,getMass)
  PYTHON_FUNC(ObjMesh,subdivide)
  PYTHON_FUNC(ObjMesh,subdivideSingle3D)
  PYTHON_FUNC(ObjMesh,subdivideSingle2D)
  PYTHON_FUNC(ObjMesh,marchingCube)
  PYTHON_FUNC(ObjMesh,makeUnique)
  PYTHON_OFUNC(ObjMesh,makeUniform,(void(ObjMesh::*)()))
  PYTHON_OFUNC(ObjMesh,makeUniform,(void(ObjMesh::*)(int,int,int,int)))
  PYTHON_FUNC(ObjMesh,smooth)
  PYTHON_FUNC(ObjMesh,insideOut)
  PYTHON_FUNC(ObjMesh,getE)
  PYTHON_FUNC(ObjMesh,buildEdge)
  PYTHON_FUNC(ObjMesh,buildKRingV)
  PYTHON_FUNC(ObjMesh,buildKRing)
  PYTHON_FUNC(ObjMesh,findInsertV)
  PYTHON_FUNC(ObjMesh,cutOpen)
  PYTHON_OFUNC(ObjMesh,getDim,(const int&(ObjMesh::*)()const))
  PYTHON_OFUNC(ObjMesh,getDim,(int&(ObjMesh::*)()));
  //-----------------------------------------------------------------Joint
  py::enum_<Joint::JOINT_TYPE>(m,"JOINT_TYPE")
  .value("TRANS_3D",Joint::TRANS_3D)
  .value("TRANS_2D",Joint::TRANS_2D)
  .value("TRANS_1D",Joint::TRANS_1D)
  .value("ROT_3D_XYZ",Joint::ROT_3D_XYZ)
  .value("ROT_3D_EXP",Joint::ROT_3D_EXP)
  .value("BALL_JOINT",Joint::BALL_JOINT)
  .value("HINGE_JOINT",Joint::HINGE_JOINT)
  .value("FIX_JOINT",Joint::FIX_JOINT)
  .value("NR_JOINT_TYPE",Joint::NR_JOINT_TYPE)
  .export_values();
  py::enum_<Joint::GEOM_TYPE>(m,"GEOM_TYPE")
  .value("MESH",Joint::MESH)
  .value("SPHERE_GEOM",Joint::SPHERE_GEOM)
  .value("SPHERE_SELF",Joint::SPHERE_SELF)
  .value("AXIS_INFO",Joint::AXIS_INFO)
  .export_values();
  py::class_<Joint>(m,"Joint")
  .def(py::init<>())
  PYTHON_FUNC(Joint,writeVTK)
  PYTHON_FUNC(Joint,getMesh)
  PYTHON_FUNC(Joint,setType)
  PYTHON_FUNC(Joint,typeToString)
  PYTHON_FUNC(Joint,assemble)
  PYTHON_OFUNC(Joint,transformMesh,(void(Joint::*)(const Mat3X4d&)))
  PYTHON_OFUNC(Joint,transformMesh,(void(Joint::*)(const Mat3d&,const Vec3d&)))
  PYTHON_OFUNC(Joint,transformMesh,(void(Joint::*)(const Vec3d&)))
  PYTHON_FUNC(Joint,CBegEnd)
  PYTHON_FUNC(Joint,RBegEnd)
  PYTHON_FUNC(Joint,nrDOF)
  PYTHON_FUNC(Joint,nrDDT)
  PYTHON_OFUNC(Joint,getMassC,(Mat6d(Joint::*)(const Mat3d&)const))
  PYTHON_OFUNC(Joint,getMassC,(Mat6d(Joint::*)()const))
  PYTHON_OFUNC(Joint,getMass,(Mat6d(Joint::*)(const Mat3d&)const))
  PYTHON_OFUNC(Joint,getMass,(Mat6d(Joint::*)()const))
  PYTHON_FUNC(Joint,getC)
  .def(STRINGIFY_OMP(isRotational),&Joint::isRotational)
  PYTHON_FUNC(Joint,isRoot)
  PYTHON_OFUNC(Joint,getSpheres,(Mat3Xd(Joint::*)(const Mat3X4d&)const))
  PYTHON_OFUNC(Joint,getSpheres,(void(Joint::*)(const Mat3X4d&,Joint::Mat3XTM)const))
  PYTHON_FUNC(Joint,avgRadGeomColl)
  PYTHON_FUNC(Joint,avgRadSelfColl)
  PYTHON_FUNC(Joint,loopAllJointTypes)
  PYTHON_MEMBER_READ(Joint,_children)
  PYTHON_MEMBER_READ(Joint,_parent)
  PYTHON_MEMBER_READ(Joint,_depth)
  PYTHON_MEMBER_READ(Joint,_typeJoint)
  PYTHON_MEMBER_READ(Joint,_mimic)
  PYTHON_MEMBER_READ(Joint,_offDOF)
  PYTHON_MEMBER_READ(Joint,_offDDT)
  PYTHON_MEMBER_READ(Joint,_limits)
  PYTHON_MEMBER_READ(Joint,_control)
  PYTHON_MEMBER_READ(Joint,_damping)
  PYTHON_MEMBER_READ(Joint,_trans)
  PYTHON_MEMBER_READ(Joint,_mult)
  PYTHON_MEMBER_READ(Joint,_offset)
  PYTHON_MEMBER_READ(Joint,_M)
  PYTHON_MEMBER_READ(Joint,_MC)
  PYTHON_MEMBER_READ(Joint,_MCCT)
  PYTHON_MEMBER_READ(Joint,_name)
  PYTHON_MEMBER_READ(Joint,_spheres)
  PYTHON_MEMBER_READ(Joint,_radGeomColl)
  PYTHON_MEMBER_READ(Joint,_radSelfColl)
  PYTHON_MEMBER_READ(Joint,_color);
  //-----------------------------------------------------------------ArticulatedBody
  py::class_<ArticulatedBody>(m,"ArticulatedBody")
  .def(py::init<>())
  .def(py::init<const tinyxml2::XMLElement&>())
  PYTHON_FUNC(ArticulatedBody,updateGeom)
  PYTHON_OFUNC(ArticulatedBody,getGeom,(const StaticGeom&(ArticulatedBody::*)()const))
  PYTHON_OFUNC(ArticulatedBody,getGeom,(StaticGeom&(ArticulatedBody::*)()))
  PYTHON_OFUNC(ArticulatedBody,getGeomEnv,(const StaticGeom&(ArticulatedBody::*)()const))
  PYTHON_OFUNC(ArticulatedBody,getGeomEnv,(StaticGeom&(ArticulatedBody::*)()))
  PYTHON_FUNC(ArticulatedBody,randomize)
  PYTHON_FUNC(ArticulatedBody,writeMesh)
  PYTHON_OFUNC(ArticulatedBody,writeVTK,(void(ArticulatedBody::*)(const Mat3Xd&,VTKWriter<scalar>&,Joint::GEOM_TYPE,const std::set<sizeType>*,const std::vector<std::vector<unsigned char> >*)const))
  PYTHON_OFUNC(ArticulatedBody,writeVTK,(void(ArticulatedBody::*)(const Mat3Xd&,const std::string&,Joint::GEOM_TYPE,const std::set<sizeType>*,const std::vector<std::vector<unsigned char> >*)const))
  PYTHON_FUNC(ArticulatedBody,writePov)
  PYTHON_FUNC(ArticulatedBody,children)
  PYTHON_FUNC(ArticulatedBody,commonRoot)
  PYTHON_FUNC(ArticulatedBody,hasJoint)
  PYTHON_FUNC(ArticulatedBody,isLeaf)
  PYTHON_FUNC(ArticulatedBody,control)
  PYTHON_FUNC(ArticulatedBody,damping)
  PYTHON_FUNC(ArticulatedBody,coefLimit)
  PYTHON_OFUNC(ArticulatedBody,lowerLimit,(Cold(ArticulatedBody::*)()const))
  PYTHON_OFUNC(ArticulatedBody,upperLimit,(Cold(ArticulatedBody::*)()const))
  PYTHON_OFUNC(ArticulatedBody,lowerLimit,(Cold(ArticulatedBody::*)(scalarD)const))
  PYTHON_OFUNC(ArticulatedBody,upperLimit,(Cold(ArticulatedBody::*)(scalarD)const))
  PYTHON_FUNC(ArticulatedBody,clampLimit)
  PYTHON_FUNC(ArticulatedBody,randomPose)
  PYTHON_FUNC(ArticulatedBody,getT)
  PYTHON_FUNC(ArticulatedBody,mimic)
  PYTHON_FUNC(ArticulatedBody,movable)
  PYTHON_OFUNC(ArticulatedBody,joint,(const Joint&(ArticulatedBody::*)(sizeType) const))
  PYTHON_OFUNC(ArticulatedBody,joint,(Joint&(ArticulatedBody::*)(sizeType)))
  PYTHON_FUNC(ArticulatedBody,jointId)
  PYTHON_FUNC(ArticulatedBody,rootJointId)
  PYTHON_FUNC(ArticulatedBody,depth)
  PYTHON_FUNC(ArticulatedBody,nrDOF)
  PYTHON_FUNC(ArticulatedBody,nrDDT)
  PYTHON_FUNC(ArticulatedBody,nrJ)
  PYTHON_FUNC(ArticulatedBody,fillChildren)
  PYTHON_FUNC(ArticulatedBody,debugBase)
  PYTHON_FUNC(ArticulatedBody,addBase)
  PYTHON_FUNC(ArticulatedBody,simplify)
  PYTHON_FUNC(ArticulatedBody,eliminateJoint)
  PYTHON_FUNC(ArticulatedBody,scaleMass)
  PYTHON_FUNC(ArticulatedBody,totalMass);
  //-----------------------------------------------------------------ArticulatedLoader
  py::class_<ArticulatedLoader>(m,"ArticulatedLoader")
  PYTHON_SFUNC(ArticulatedLoader,readURDF)
  PYTHON_SFUNC(ArticulatedLoader,compareVisualMeshURDF)
  PYTHON_SFUNC(ArticulatedLoader,visualizeJointLimitVTK)
  //build-in bodies
  PYTHON_OSFUNC(ArticulatedLoader,createBird,(ArticulatedBody(*)(sizeType,bool)))
  PYTHON_OSFUNC(ArticulatedLoader,createBipedal,(ArticulatedBody(*)(sizeType,bool)))
  PYTHON_OSFUNC(ArticulatedLoader,createChain,(ArticulatedBody(*)(sizeType,scalar,sizeType,sizeType)))
  PYTHON_OSFUNC(ArticulatedLoader,createSpider,(ArticulatedBody(*)(sizeType,scalar,scalar)))
  PYTHON_OSFUNC(ArticulatedLoader,createArm,(ArticulatedBody(*)(scalar,scalar)));
  //-----------------------------------------------------------------ExternalWrench
#define EXTERNAL_WRENCH(T,SUFFIX)   \
py::class_<ExternalWrench<T>>(m,("ExternalWrench"+std::string(SUFFIX)).c_str()) \
.def(py::init<>())  \
PYTHON_FIELD_READWRITE(ExternalWrench<T>,_DBDq) \
PYTHON_FIELD_READWRITE(ExternalWrench<T>,_DBDX) \
PYTHON_FIELD_READWRITE(ExternalWrench<T>,_jointId)  \
PYTHON_FIELD_READWRITE(ExternalWrench<T>,_B)    \
PYTHON_FIELD_READWRITE(ExternalWrench<T>,_w);
  EXTERNAL_WRENCH(double,"d")
  EXTERNAL_WRENCH(__float128,"q")
  EXTERNAL_WRENCH(mpfr::mpreal,"m")
  //-----------------------------------------------------------------MDP
#define MDP(T,SUFFIX)   \
py::class_<MDP<T>>(m,("MDP"+std::string(SUFFIX)).c_str())   \
.def(py::init<const ArticulatedBody&,Options&>())   \
.def("solveMDPQPF",[](MDP<T>& mdp,MDP<T>::DMat* DdqHatDq,MDP<T>::DMat* DdqHatDdq,bool profileMDPError) {   \
  return mdp.solveMDPQPF(mdp.mapM(DdqHatDq),mdp.mapM(DdqHatDdq),profileMDPError);  \
})  \
.def("solveMDPQPI",[](MDP<T>& mdp,MDP<T>::DMat* DdqHatDq,MDP<T>::DMat* DdqHatDdq,bool profileMDPError) {   \
  return mdp.solveMDPQPI(mdp.mapM(DdqHatDq),mdp.mapM(DdqHatDdq),profileMDPError);  \
})  \
.def("solveMDPLPF",[](MDP<T>& mdp,MDP<T>::DMat* DddqDq,MDP<T>::DMat* DddqDdq,bool profileMDPError) {   \
  return mdp.solveMDPQPF(mdp.mapM(DddqDq),mdp.mapM(DddqDdq),profileMDPError);  \
})  \
.def("solveMDPLPI",[](MDP<T>& mdp,MDP<T>::DMat* DddqDq,MDP<T>::DMat* DddqDdq,bool profileMDPError) {   \
  return mdp.solveMDPQPI(mdp.mapM(DddqDq),mdp.mapM(DddqDdq),profileMDPError);  \
})  \
PYTHON_FUNC(MDP<T>,reset)   \
PYTHON_FUNC(MDP<T>,randomize)   \
PYTHON_FUNC(MDP<T>,readAndTestProb) \
PYTHON_FUNC(MDP<T>,assembleHcQP)  \
PYTHON_FUNC(MDP<T>,assembleHcLP)  \
PYTHON_FUNC(MDP<T>,B)   \
PYTHON_FUNC(MDP<T>,Bf)   \
PYTHON_FUNC(MDP<T>,BT)  \
PYTHON_FUNC(MDP<T>,DBDq)    \
PYTHON_FUNC(MDP<T>,DBDqT)   \
PYTHON_FUNC(MDP<T>,DXDq)    \
PYTHON_FUNC(MDP<T>,kineticEnergyInfo)   \
PYTHON_FUNC(MDP<T>,kineticEnergyAssembled)  \
PYTHON_SFUNC(MDP<T>,debugWrenchConstructor) \
PYTHON_SFUNC(MDP<T>,debug)  \
PYTHON_FIELD_READWRITE(MDP<T>,_dqMDP)   \
PYTHON_FIELD_READWRITE(MDP<T>,_lastW)   \
PYTHON_FIELD_READWRITE(MDP<T>,_externalWrench);
  MDP(double,"d")
  MDP(__float128,"q")
  MDP(mpfr::mpreal,"m")
  //-----------------------------------------------------------------C0EnvWrenchConstructor
#define C0_ENV_WRENCH_CONSTRUCTOR(T,SUFFIX)   \
py::class_<C0EnvWrenchConstructor<T>,std::shared_ptr<C0EnvWrenchConstructor<T>>>(m,("C0EnvWrenchConstructor"+std::string(SUFFIX)).c_str())    \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,scalarD,scalarD,scalarD,sizeType,const C0EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T>()) \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,std::function<scalarD(scalarD,scalarD)>,sizeType,const C0EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T>()) \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,scalarD,sizeType,scalarD,const C0EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T>()) \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,scalarD,const C0EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T>()) \
.def(py::init<const ArticulatedBody&,const Vec4d&,const C0EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T>()) \
.def(py::init<const ArticulatedBody&,const std::string&,bool,scalarD,const C0EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T>()) \
PYTHON_FUNC(C0EnvWrenchConstructor<T>,writeEndEffectorVTK)    \
PYTHON_FUNC(C0EnvWrenchConstructor<T>,writeVTK)   \
PYTHON_FIELD_READWRITE(C0EnvWrenchConstructor<T>,_externalForces) \
PYTHON_MEMBER_READ(C0EnvWrenchConstructor<T>,_frm)   \
PYTHON_MEMBER_READ(C0EnvWrenchConstructor<T>,_B);
  C0_ENV_WRENCH_CONSTRUCTOR(double,"d")
  C0_ENV_WRENCH_CONSTRUCTOR(__float128,"q")
  C0_ENV_WRENCH_CONSTRUCTOR(mpfr::mpreal,"m")
  //-----------------------------------------------------------------C2EnvWrenchConstructor
#define C2_ENV_WRENCH_CONSTRUCTOR(T,SUFFIX)   \
py::class_<C2EnvWrenchConstructor<T>,C0EnvWrenchConstructor<T>,std::shared_ptr<C2EnvWrenchConstructor<T>>>(m,("C2EnvWrenchConstructor"+std::string(SUFFIX)).c_str())    \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,scalarD,scalarD,scalarD,sizeType,const C2EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T,T>()) \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,std::function<scalarD(scalarD,scalarD)>,sizeType,const C2EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T,T>()) \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,scalarD,sizeType,scalarD,const C2EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T,T>()) \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,scalarD,const C2EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T,T>()) \
.def(py::init<const ArticulatedBody&,const Vec4d&,const C2EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T,T>()) \
.def(py::init<const ArticulatedBody&,const std::string&,bool,scalarD,const C2EnvWrenchConstructor<T>::Vec3T&,sizeType,T,T,T>()) \
PYTHON_MEMBER_READ(C2EnvWrenchConstructor<T>,_crossB)   \
PYTHON_MEMBER_READ(C2EnvWrenchConstructor<T>,_coef);
  C2_ENV_WRENCH_CONSTRUCTOR(double,"d")
  C2_ENV_WRENCH_CONSTRUCTOR(__float128,"q")
  C2_ENV_WRENCH_CONSTRUCTOR(mpfr::mpreal,"m")
  //-----------------------------------------------------------------C2GranularWrenchConstructor
#define C2_GRANULAR_WRENCH_CONSTRUCTOR(T,SUFFIX)   \
py::class_<C2GranularWrenchConstructor<T>,std::shared_ptr<C2GranularWrenchConstructor<T>>>(m,("C2GranularWrenchConstructor"+std::string(SUFFIX)).c_str())    \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,scalarD,scalarD,scalarD,sizeType,const std::string&,const std::string&>()) \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,std::function<scalarD(scalarD,scalarD)>,sizeType,const std::string&,const std::string&>()) \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,scalarD,sizeType,scalarD,const std::string&,const std::string&>()) \
.def(py::init<const ArticulatedBody&,scalarD,scalarD,scalarD,const std::string&,const std::string&>()) \
.def(py::init<const ArticulatedBody&,const Vec4d&,const std::string&,const std::string&>())   \
.def(py::init<const ArticulatedBody&,const std::string&,bool,scalarD,const std::string&,const std::string&>())   \
PYTHON_FUNC(C2GranularWrenchConstructor<T>,writeEndEffectorVTK)    \
PYTHON_FUNC(C2GranularWrenchConstructor<T>,writeVTK);
  C2_GRANULAR_WRENCH_CONSTRUCTOR(double,"d")
  C2_GRANULAR_WRENCH_CONSTRUCTOR(__float128,"q")
  C2_GRANULAR_WRENCH_CONSTRUCTOR(mpfr::mpreal,"m")
  //-----------------------------------------------------------------PDTarget
  py::class_<PDTarget,std::shared_ptr<PDTarget>>(m,"PDTarget")
      .def(py::init<const PDTarget::Vec&,const PDTarget::Vec&,const std::vector<std::tuple<scalarD,PDTarget::Vec,PDTarget::Vec>>&>())
      .def(py::init<const PDTarget::Vec&,const PDTarget::Vec&,const PDTarget::Vec&>())
      PYTHON_FUNC(PDTarget,writeVTKSeq)
      PYTHON_FUNC(PDTarget,reset)
      PYTHON_FUNC(PDTarget,advance)
      PYTHON_FUNC(PDTarget,s)
      PYTHON_FUNC(PDTarget,setS)
      PYTHON_FUNC(PDTarget,splineInterp)
      PYTHON_FIELD_READWRITE(PDTarget,_PCoef)
      PYTHON_FIELD_READWRITE(PDTarget,_DCoef);
  //-----------------------------------------------------------------MDPSimulator
  py::enum_<MDP_SIMULATOR_MODE>(m,"MDP_SIMULATOR_MODE")
  .value("INVERSE_I",MDP_SIMULATOR_MODE::INVERSE_I)
  .value("FORWARD_I",MDP_SIMULATOR_MODE::FORWARD_I)
  .value("INVERSE_LF",MDP_SIMULATOR_MODE::INVERSE_LF)
  .value("FORWARD_LF",MDP_SIMULATOR_MODE::FORWARD_LF)
  .value("INVERSE_BACKWARD_RK1F",MDP_SIMULATOR_MODE::INVERSE_BACKWARD_RK1F)
  .value("BACKWARD_RK1F",MDP_SIMULATOR_MODE::BACKWARD_RK1F)
  .value("FORWARD_RK1F",MDP_SIMULATOR_MODE::FORWARD_RK1F)
  .value("FORWARD_RK1I",MDP_SIMULATOR_MODE::FORWARD_RK1I)
  .value("FORWARD_RK2I",MDP_SIMULATOR_MODE::FORWARD_RK2I)
  .value("FORWARD_RK4I",MDP_SIMULATOR_MODE::FORWARD_RK4I)
  .value("NMDP_PGM",MDP_SIMULATOR_MODE::NMDP_PGM)
  .value("NMDP_ZOPGM",MDP_SIMULATOR_MODE::NMDP_ZOPGM)
  .export_values();
#define MDP_SIMULATOR(T,SUFFIX)   \
py::class_<MDPSimulator<T>,MDP<T>>(m,("MDPSimulator"+std::string(SUFFIX)).c_str())    \
.def(py::init<const ArticulatedBody&,Options&,const MDPSimulator<T>::Vec3T&,MDP_SIMULATOR_MODE>())    \
PYTHON_FUNC(MDPSimulator<T>,setWrenchConstructor)   \
.def("setC0EnvWrenchConstructor",[](MDPSimulator<T>& v,std::shared_ptr<C0EnvWrenchConstructor<T>> f){v.setWrenchConstructor(f);})  \
.def("setC2EnvWrenchConstructor",[](MDPSimulator<T>& v,std::shared_ptr<C2EnvWrenchConstructor<T>> f){v.setWrenchConstructor(f);})  \
.def("setC2GranularWrenchConstructor",[](MDPSimulator<T>& v,std::shared_ptr<C2GranularWrenchConstructor<T>> f){v.setWrenchConstructor(f);})  \
.def("getTerrainMesh",[](const MDPSimulator<T>& v) {    \
  MDPSimulator<T>::WrenchConstructor w=v.getWrenchConstructor();    \
  std::shared_ptr<EnvWrenchConstructor<T>> env=std::dynamic_pointer_cast<EnvWrenchConstructor<T>>(w);   \
  if(env)   \
    return env->_env->getMesh();  \
  else return ObjMesh();    \
})  \
PYTHON_FUNC(MDPSimulator<T>,setPDTarget)    \
PYTHON_FUNC(MDPSimulator<T>,writeVTKSeq)    \
PYTHON_FUNC(MDPSimulator<T>,writeVTK)   \
PYTHON_FUNC(MDPSimulator<T>,setState)   \
.def("step",[](MDPSimulator<T>& sim,T dt,const MDPSimulator<T>::Vec* tau0,const MDPSimulator<T>::Vec& qdq,bool Dqdq,bool Dtau){  \
  MDPSimulator<T>::DMat Dqdqs,Dtaus;    \
  MDPSimulator<T>::Vec o=sim.step(dt,tau0,qdq,Dqdq?&Dqdqs:NULL,Dtau?&Dtaus:NULL);    \
  return py::make_tuple(o,Dqdqs,Dtaus); \
})  \
.def("stepBatched",[](MDPSimulator<T>& sim,T dt,const MDPSimulator<T>::Vecss* tau0,const MDPSimulator<T>::Vecss& qdq,bool Dqdq,bool Dtau){   \
  MDPSimulator<T>::DMatss Dqdqs,Dtaus;    \
  if(Dqdq)Dqdqs.resize(qdq.size()); \
  if(Dtau)Dtaus.resize(qdq.size()); \
  MDPSimulator<T>::Vecss o=sim.stepBatched(dt,tau0,qdq,Dqdq?&Dqdqs:NULL,Dtau?&Dtaus:NULL);  \
  return py::make_tuple(o,Dqdqs,Dtaus);   \
})  \
.def("getJointTransform",[](MDPSimulator<T>& sim,sizeType id)->MDPSimulator<T>::Mat3XT {    \
  if(id<0 || id>=sim._body.nrJ())   \
    return MDPSimulator<T>::Mat3XT::Zero(3,0);   \
  else return TRANSI(sim._info.getTrans(),id);  \
})  \
.def("getJointTransform",[](MDPSimulator<T>& sim,const std::string& name)->MDPSimulator<T>::Mat3XT {    \
  sizeType id=-1;   \
  for(sizeType i=0;i<sim._body.nrJ();i++)    \
    if(sim._body.joint(i)._name==name) {    \
      id=i; \
      break;    \
    }   \
  if(id<0 || id>=sim._body.nrJ())   \
    return MDPSimulator<T>::Mat3XT::Zero(3,0);   \
  else return TRANSI(sim._info.getTrans(),id);  \
})  \
.def("getJointForce",[](MDPSimulator<T>& sim) { \
  std::vector<std::pair<sizeType,MDPSimulator<T>::Vec3T>> ret;    \
  for(sizeType i=0;i<(sizeType)sim._externalWrench.size();i++) {  \
    MDPSimulator<T>::Vec3T f;   \
    if(sim._externalWrench[i]._w.size()==0) \
      f.setZero();  \
    else f=(sim._externalWrench[i]._B*sim._externalWrench[i]._w).template segment<3>(3);    \
    ret.push_back(std::make_pair(sim._externalWrench[i]._jointId,f));    \
  } \
  return ret;   \
})  \
PYTHON_FUNC(MDPSimulator<T>,debug)  \
PYTHON_FUNC(MDPSimulator<T>,setMode)\
PYTHON_FUNC(MDPSimulator<T>,getMode)\
PYTHON_FUNC(MDPSimulator<T>,begLog) \
PYTHON_FUNC(MDPSimulator<T>,endLog);
  MDP_SIMULATOR(double,"d")
  MDP_SIMULATOR(__float128,"q")
  MDP_SIMULATOR(mpfr::mpreal,"m")
  //-----------------------------------------------------------------EndEffectorBounds
  py::class_<EndEffectorBounds>(m,"EndEffectorBounds")
  .def(py::init<>())
  .def(py::init<sizeType>())
  .def(py::init<sizeType,const Vec3d&,scalarD>())
  PYTHON_FUNC(EndEffectorBounds,exists)
  PYTHON_FUNC(EndEffectorBounds,isBetterBB)
  PYTHON_FUNC(EndEffectorBounds,sameOrthant)
  PYTHON_FUNC(EndEffectorBounds,closeFar)
  PYTHON_FUNC(EndEffectorBounds,nrDOF)
  PYTHON_FUNC(EndEffectorBounds,getDOF)
  PYTHON_FUNC(EndEffectorBounds,globalPosAll)
  PYTHON_FUNC(EndEffectorBounds,globalPos)
  PYTHON_FUNC(EndEffectorBounds,randomPos)
  PYTHON_FUNC(EndEffectorBounds,jointId)
  .def("detectEndEffector",[](EndEffectorBounds& ee,const ArticulatedBody& body,const Vec3d& zRange) {
    SimplifiedDynamics::detectEndEffector(body,ee._JID.back(),ee._localPos,ee._phi0,zRange);
  })
  PYTHON_FIELD_READWRITE(EndEffectorBounds,_JID)
  PYTHON_FIELD_READWRITE(EndEffectorBounds,_bb)
  PYTHON_FIELD_READWRITE(EndEffectorBounds,_localPos)
  PYTHON_FIELD_READWRITE(EndEffectorBounds,_phi0)
  PYTHON_FIELD_READWRITE(EndEffectorBounds,_res);
  //-----------------------------------------------------------------Options
  py::class_<Options>(m,"Options")
  .def(py::init<>())
  .def("setOption",[](Options& op,const MDPSimulator<double>&,const std::string& name,double val) {
    op.setOptions<MultiPrecisionLQP<double>,double>(name,val);
  })
  .def("setOption",[](Options& op,const MDPSimulator<__float128>&,const std::string& name,__float128 val) {
    op.setOptions<MultiPrecisionLQP<__float128>,__float128>(name,val);
  })
  .def("setOption",[](Options& op,const MDPSimulator<mpfr::mpreal>&,const std::string& name,mpfr::mpreal val) {
    op.setOptions<MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal>(name,val);
  })
  .def("setOption",[](Options& op,const MDPSimulator<double>&,const std::string& name,bool val) {
    op.setOptions<MultiPrecisionLQP<double>,bool>(name,val);
  })
  .def("setOption",[](Options& op,const MDPSimulator<__float128>&,const std::string& name,bool val) {
    op.setOptions<MultiPrecisionLQP<__float128>,bool>(name,val);
  })
  .def("setOption",[](Options& op,const MDPSimulator<mpfr::mpreal>&,const std::string& name,bool val) {
    op.setOptions<MultiPrecisionLQP<mpfr::mpreal>,bool>(name,val);
  })
  .def("setOption",[](Options& op,const MDPSimulator<double>&,const std::string& name,int val) {
    op.setOptions<MultiPrecisionLQP<double>,int>(name,val);
  })
  .def("setOption",[](Options& op,const MDPSimulator<__float128>&,const std::string& name,int val) {
    op.setOptions<MultiPrecisionLQP<__float128>,int>(name,val);
  })
  .def("setOption",[](Options& op,const MDPSimulator<mpfr::mpreal>&,const std::string& name,int val) {
    op.setOptions<MultiPrecisionLQP<mpfr::mpreal>,int>(name,val);
  })
  .def("printVars",[](const Options& v) {
    v.print();
  });
}
