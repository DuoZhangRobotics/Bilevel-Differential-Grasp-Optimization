#include "graspit/EGPlanner/energy/contactEnergy.h"
#include "graspit/robot.h"
#include "graspit/grasp.h"
#include "graspit/debug.h"
#include "graspit/world.h"
#include "graspit/contact/virtualContact.h"

double ContactEnergy::energy() const
{
  /*
  this is if we don't want to use pre-specified contacts, but rather the closest points between each link
  and the object all iterations. Not necessarily needed for contact energy, but useful for GQ energy
  */

  if (mContactType == CONTACT_LIVE && mType != "AUTO_GRASP_QUALITY_ENERGY" && mType != "STRICT_AUTO_GRASP_ENERGY")
  {
    mHand->getWorld()->findVirtualContacts(mHand, mObject);
    DBGP("Live contacts computation");
  }

  mHand->getGrasp()->collectVirtualContacts();

  //DBGP("Contact energy computation")
  //average error per contact
  VirtualContact *contact;
  vec3 p, n, cn;
  double totalError = 0;
  for (int i = 0; i < mHand->getGrasp()->getNumContacts(); i++)
  {
    contact = (VirtualContact *)mHand->getGrasp()->getContact(i);
    contact->getObjectDistanceAndNormal(mObject, &p, NULL);
    double dist = p.norm();

    //this should never happen anymore since we're never inside the object
    //if ( (-1.0 * p) % n < 0) dist = -dist;

    //BEST WORKING VERSION, strangely enough
    totalError += fabs(dist);

    //let's try this some more
    //totalError += distanceFunction(dist);
    //cn = -1.0 * contact->getWorldNormal();

    //new version
    cn = contact->getWorldNormal();
    n = p.normalized();
    double d = 1 - cn.dot(n);
    totalError += d * 100.0 / 2.0;
  }

  totalError /= mHand->getGrasp()->getNumContacts();

  //DBGP("Contact energy: " << totalError);
  return totalError;
}

//=======================================================================================================================
// #include <Python.h>
// #include <memory>
// #include "graspit/EGPlanner/energy/contactEnergy.h"
// #include "graspit/robot.h"
// #include "graspit/grasp.h"
// #include "graspit/debug.h"
// #include "graspit/world.h"
// #include "graspit/contact/virtualContact.h"
// #include "graspit/quality/quality.h"
// #include <iostream>
// #include <string.h>

/**
 * @brief This formulation combines virtual contact energy with autograsp energy. In addition, it computes EBM-based energy to access human-like grasping. 
    Virtual contact energy is used to "guide" initial stages of the search and to see if we should even bother computing autograsp quality. 
    Autograsp is a couple of orders of magnitude higher and so should work very well with later stages of the sim ann search.
    EBM-based energy is used to guide human-like grasping generation. This energy is computed by EBM learning model.
    Note that, EBM-based energy is implementation based on python, so we need to write an interface to run python script.
 * @author Jian Liu
 * 
 */

// PyObject *pModule = NULL;
// PyObject *pFunc = NULL;
// PyObject *pReturn = NULL;
// PyObject *list_dof = NULL;
// PyObject *pArgs = NULL;

// double
// ContactEnergy::energy() const
// {
//   //first compute regular contact energy; also count how many links are "close" to the object
//   if (mContactType == CONTACT_LIVE && mType != "AUTO_GRASP_QUALITY_ENERGY" && mType != "STRICT_AUTO_GRASP_ENERGY")
//   {
//     mHand->getWorld()->findVirtualContacts(mHand, mObject);
//     DBGP("Live contacts computation");
//   }

//   mHand->getGrasp()->collectVirtualContacts();

//   //DBGP("Contact energy computation")
//   //average error per contact
//   VirtualContact *contact;
//   vec3 p, n, cn;
//   double virtualError = 0.0;
//   int closeContacts = 0;
//   double ebmQuality = 0.0;
//   int numDof = mHand->getNumDOF();
//   for (int i = 0; i < mHand->getGrasp()->getNumContacts(); i++)
//   {
//     contact = (VirtualContact *)mHand->getGrasp()->getContact(i);
//     contact->getObjectDistanceAndNormal(mObject, &p, NULL);
//     double dist = p.norm();

//     //BEST WORKING VERSION, strangely enough
//     virtualError += fabs(dist);

//     //new version
//     cn = contact->getWorldNormal();
//     n = p.normalized();
//     double d = 1 - cn.dot(n);
//     virtualError += d * 100.0 / 2.0;

//     if (fabs(dist) < 20 && d < 0.3)
//     {
//       closeContacts++;
//     }
//   }

//   virtualError /= mHand->getGrasp()->getNumContacts();

//   double *dofVals = new double[numDof];
//   mHand->getDOFVals(dofVals);

//   std::vector<double> DOFs(numDof);
//   for (int i = 0; i < numDof; i++)
//   {
//     DOFs[i] = dofVals[i];
//   }
//   delete[] dofVals;
//   ebmQuality = this->method1(DOFs, numDof, mHand->getEBMPath());

//   //Test
//   // std::cout<<mHand->getNumDOF()<<std::endl;
//   // for(int i=0;i<mHand->getNumDOF();i++){
//   //   std::cout<<dofVals[i]<<std::endl;
//   //   }
//   // std::cout<<mHand->getEBMPath()<<std::endl;
//   // if (this->Initialize(mHand->getEBMPath()) < 0)
//   // {
//   //   std::cout<< "Initialize is failed"<<std::endl;
//   // }
//   // ebmQuality = this->ebm_pythonInterface(DOFs,numDof);
//   // //virtualError = virtualError + ebmQuality * 1.0e3;
//   // this->Uninitialize();

//   //if more than 2 links are "close" go ahead and compute the true quality
//   double volQuality = 0.0, epsQuality = 0.0;
//   if (closeContacts >= 2)
//   {
//     mHand->autoGrasp(false, 1.0);
//     //now collect the true contacts;
//     mHand->getGrasp()->collectContacts();
//     if (mHand->getGrasp()->getNumContacts() >= 4)
//     {
//       mHand->getGrasp()->updateWrenchSpaces();
//       volQuality = mVolQual->evaluate();
//       epsQuality = mEpsQual->evaluate();
//       if (epsQuality < 0)
//       {
//         epsQuality = 0;
//       } //QM returns -1 for non-FC grasps
//     }

//     DBGP("Virtual error " << virtualError << " and " << closeContacts << " close contacts.");
//     DBGP("Volume quality: " << volQuality << " Epsilon quality: " << epsQuality);
//     DBGP("Human-like quality: " << ebmQuality);
//   }
//   std::cout << "contact energy: " << virtualError << std::endl;
//   std::cout << "Volume energy: " << volQuality * 1.0e3 << std::endl;
//   std::cout << "ebm energy: " << ebmQuality * 1.0e2 << std::endl;

//   return virtualError;
//   // if (volQuality == 0) { q = virtualError + ebmQuality * 1.0e3; }
//   // else { q = virtualError - volQuality * 1.0e3 + ebmQuality * 1.0e3; }
//   // if (volQuality || epsQuality) {DBGP("Final quality: " << q);}

//   // //DBGP("Final value: " << q << std::endl);
//   // return q;
// }

// double
// ContactEnergy::method1(std::vector<double> &dofVals, int numDOF, std::string modelPath) const
// {
//   PyRun_SimpleString("import sys");
//   //The file path of the ebmPythonInterface.py
//   PyRun_SimpleString("sys.path.append('/root/WorkSpace/EBM_Hand')");
//   pModule = PyImport_ImportModule("ebmPythonInterface");

//   if (!pModule)
//   {
//     std::cout << "pModle is null" << std::endl;
//   }

//   if (PyImport_ImportModule("ebmPythonInterface") == NULL || PyErr_Occurred())
//   {
//     PyErr_Print();
//   }

//   while (!pModule)
//   {
//     sleep(3.0);
//   }

//   if (!pModule)
//   {
//     std::cout << "Error: python module is null!" << std::endl;
//     return 0;
//   }
//   else
//   {
//     std::cout << "python module is successful" << std::endl;
//   }

//   pFunc = PyObject_GetAttrString(pModule, "dof_ebm");
//   if (!pFunc)
//   {
//     std::cout << "Error: python pFunc_dof_ebm is null!" << std::endl;
//     return 0;
//   }
//   else
//   {
//     std::cout << "python pFunc_dof_ebm is successful" << std::endl;
//   }

//   //创建参数:
//   pArgs = PyTuple_New(2); //函数调用的参数传递均是以元组的形式打包的,2表示参数个数
//   list_dof = PyList_New(0);

//   for (int i = 0; i < numDOF; i++)
//   {
//     double degVal = dofVals[i] * 57.3;
//     //std::cout<<degVal<<std::endl;
//     PyList_Append(list_dof, Py_BuildValue("d", degVal));
//   }

//   PyTuple_SetItem(pArgs, 0, list_dof);

//   //char* p = new char[100];p = strcpy(p,modelPath.c_str());
//   //std::cout<<p<<std::endl;
//   PyTuple_SetItem(pArgs, 1, Py_BuildValue("s", modelPath.c_str()));

//   //返回值
//   pReturn = PyEval_CallObject(pFunc, pArgs); //调用函数
//   //将返回值转换为double类型
//   double result = 0.0;
//   PyArg_Parse(pReturn, "d", &result); //d表示转换成double型变量
//   std::cout << "ebm value: " << result << std::endl;

//   Py_DECREF(pArgs);
//   Py_DECREF(pModule);
//   Py_DECREF(pFunc);
//   Py_DECREF(list_dof);
//   Py_DECREF(pReturn);

//   pModule = NULL;
//   pFunc = NULL;
//   pReturn = NULL;
//   pArgs = NULL;
//   list_dof = NULL;
//   return result;
// }

//=======================================================================================================================