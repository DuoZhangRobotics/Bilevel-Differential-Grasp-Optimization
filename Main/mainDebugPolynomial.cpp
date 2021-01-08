#include <Deformable/SoftFinger.h>
#include <Deformable/FEMReducedSystem.h>
#include <Deformable/FEMEnergy.h>
#include <Deformable/FEMMesh.h>
#include <Utils/Hash.h>

USE_PRJ_NAMESPACE

template <typename T>
void debug()
{
  std::unordered_set<Vec3i,Hash> css;
  for(sizeType x=0; x<2; x++)
    for(sizeType y=0; y<2; y++)
      for(sizeType z=0; z<2; z++)
        css.insert(Vec3i(x,y,z));
  writeABQVoxel(css,"testPolynomial.abq",1,true);

  tinyxml2::XMLDocument doc;
  doc.InsertEndChild(doc.NewElement("root"));
  put<std::string>(doc,"path","testPolynomial.abq");
  std::shared_ptr<FEMMesh<T>> mesh(new FEMMesh<T>);
  mesh->reset(doc);
  std::shared_ptr<FEMReducedSystem<T>> systemFEM(new FEMReducedSystem<T>(mesh));
  systemFEM->addStVKElastic(10,100);
  systemFEM->addPressure(0,1,FEMMesh<T>::Vec3T::Zero(),0,0);
  systemFEM->B().setRandom(systemFEM->proj().cols(),10);
  systemFEM->PB()=systemFEM->proj()*systemFEM->B();

  FEMPolynomialEnergy<T>(*systemFEM,[&](const FEMEnergy<T>& e) {
    return dynamic_cast<const FEMStVKElasticEnergy<T>*>(&e)!=NULL;
  }).debug(*systemFEM,[&](const FEMEnergy<T>& e) {
    return dynamic_cast<const FEMStVKElasticEnergy<T>*>(&e)!=NULL;
  });

  FEMPolynomialEnergy<T>(*systemFEM,[&](const FEMEnergy<T>& e) {
    return dynamic_cast<const FEMPressureEnergy<T>*>(&e)!=NULL || dynamic_cast<const FEMEmbeddedPressureEnergy<T>*>(&e)!=NULL;
  }).debug(*systemFEM,[&](const FEMEnergy<T>& e) {
    return dynamic_cast<const FEMPressureEnergy<T>*>(&e)!=NULL || dynamic_cast<const FEMEmbeddedPressureEnergy<T>*>(&e)!=NULL;
  });
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  debug<double>();
  return 0;
}
