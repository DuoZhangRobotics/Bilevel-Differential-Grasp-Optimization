#include <Python.h>
#include <memory>
#include <iostream>
#include <string.h>
#include <cmath>
#include "graspit/EGPlanner/energy/contactEnergy.h"
#include "graspit/EGPlanner/energy/ebmGuidedAutoGraspQualityEnergy.h"
#include "graspit/robot.h"
#include "graspit/grasp.h"
#include "graspit/debug.h"
#include "graspit/world.h"
#include "graspit/contact/virtualContact.h"
#include "graspit/quality/quality.h"

/**
 * @brief This formulation combines virtual contact energy with autograsp energy. In addition, it computes EBM-based energy to access human-like grasping. 
    Virtual contact energy is used to "guide" initial stages of the search and to see if we should even bother computing autograsp quality. 
    Autograsp is a couple of orders of magnitude higher and so should work very well with later stages of the sim ann search.
    EBM-based energy is used to guide human-like grasping generation. This energy is computed by EBM learning model.
    Note that, EBM-based energy is implementation based on python, so we need to write an interface to run python script.
 * @author Jian Liu
 * 
 */

PyObject *pModule = NULL;
PyObject *pFunc = NULL;
PyObject *pReturn = NULL;
PyObject *list_dof = NULL;
PyObject *pArgs = NULL;
int humanLike_legal = -1;

// bool EBMGuidedAutoGraspQualityEnergy::legal() const
// {
//   //full collision detection
//   //if the hand is passed as an argument, this should only check for collisions that
//   //actually involve the hand
//   bool collision_flag = mHand->getWorld()->noCollision(mHand);
//   if (humanLike_legal == 1)
//   {
//     humanLike_legal = -1; //rest humanlike_legal
//     std::cout << "The current hand is human-like." << std::endl;
//     return collision_flag;
//   }
//   else if (humanLike_legal == 0)
//   {
//     humanLike_legal = -1; ////rest humanlike_legal
//     std::cout << "The current hand is not human-like." << std::endl;
//     return collision_flag;
//   }
//   else
//   {
//     if (mType == "COMPLIANT_ENERGY" || mType == "DYNAMIC_AUTO_GRASP_ENERGY")
//     {
//       return true;
//     }
//   }
//   return collision_flag;
// }

// double
// EBMGuidedAutoGraspQualityEnergy::energy() const
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
//   double AutograspError = 0.0;
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
//   double volQuality = 0.0;
//   double epsQuality = 0.0;

//   if (closeContacts >= 3) //比较接近的手上接触点的个数越多越好，人手总过有18个接触点
//   {
//     //Evaluating human-like energy of the grasp candiate after closing the human hand
//     double *dofVals = new double[numDof];
//     mHand->getDOFVals(dofVals);

//     std::vector<double> DOFs(numDof);
//     for (int i = 0; i < numDof; i++)
//     {
//       DOFs[i] = dofVals[i];
//     }
//     delete[] dofVals;
//     ebmQuality = this->ebm_pythonInterface(DOFs, numDof, mHand->getEBMPath());

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

//       //AutograspError = -(30 * volQuality) - (100 * epsQuality);
//       AutograspError = -volQuality * 1.0e3;
//     }

//     DBGP("Virtual error " << virtualError << " and " << closeContacts << " close contacts.");
//     DBGP("Volume quality: " << volQuality << " Epsilon quality: " << epsQuality);
//     DBGP("Human-like quality: " << ebmQuality);
//   }
//   std::cout << "Virtual error : " << virtualError << std::endl;
//   std::cout << "Volume quality : " << volQuality << " Epsilon quality: " << epsQuality << std::endl;
//   std::cout << "AutograspError energy: " << AutograspError << std::endl;
//   std::cout << "ebm energy: " << ebmQuality << std::endl;

//   double q;
//   if (volQuality == 0)
//   {
//     q = virtualError;
//     return q;
//   }
//   else
//   {
//     q = virtualError + AutograspError;
//   }

//   // if (volQuality || epsQuality)
//   // {
//   //   DBGP("Final quality: " << q);
//   // }

//   if (ebmQuality < 0.0)
//   {
//     // if( q < -1.0){
//     //   humanLike_legal = 0;//grasp hand is human-like
//     // }else{
//     //   humanLike_legal = 1;//grasp hand is human-like
//     // }
//     humanLike_legal = 1; //grasp hand is human-like
//     return q;
//   }
//   else
//   {
//     //research eigengrasp again.
//     humanLike_legal = 0; //grasp hand is not human-like
//     return q + ebmQuality * 200.0;
//   }
//   // double q;
//   // if (ebmQuality >= 0.5)
//   // {
//   //   if (volQuality == 0)
//   //   {
//   //     //q = virtualError + ebmQuality * 500.0;
//   //     q = virtualError * (1.0 + ebmQuality);
//   //   }
//   //   else
//   //   {
//   //     q = (virtualError + AutograspError) * (1.0 + ebmQuality);
//   //   }
//   // }
//   // else if (ebmQuality >= 0 && ebmQuality < 0.5)
//   // {
//   //   if (volQuality == 0)
//   //   {
//   //     q = virtualError * (1.0 + ebmQuality);
//   //   }
//   //   else
//   //   {
//   //     q = (virtualError + AutograspError) * (1.0 + ebmQuality);
//   //   }
//   // }
//   // else
//   // {
//   //   if (volQuality == 0)
//   //   {
//   //     if (virtualError >= 0)
//   //     {
//   //       q = -1.0 * virtualError * ebmQuality;
//   //     }
//   //     else
//   //     {
//   //       q = (1.0 - ebmQuality) * virtualError;
//   //     }
//   //   }
//   //   else
//   //   {
//   //     if (virtualError >= 0)
//   //     {
//   //       q = -1.0 * (virtualError + AutograspError) * ebmQuality;
//   //     }
//   //     else
//   //     {
//   //       q = (1.0 - ebmQuality) * (virtualError + AutograspError);
//   //     }
//   //   }
//   // }
//   // //The solid constraint of grasp energy computation
//   // if (volQuality || epsQuality)
//   // {
//   //   DBGP("Final quality: " << q);
//   // }
//   // std::cout << "Final grasp energy value: " << q << std::endl;
//   // return q;
// }

double
EBMGuidedAutoGraspQualityEnergy::energy() const
{
  int closeContacts = 0;
  double ContactError = contactEnergy(closeContacts);
  double AutograspError = approachAutograspQualityEnergy();
  double HumanLikeError = ebmEnergy();
  std::cout<<"Close-contacts: "<<closeContacts<<std::endl;

  if (closeContacts >= 2) //比较接近的手上接触点的个数越多越好，人手总过有18个接触点
  {
    if (AutograspError != 0.0)
    {
      if(HumanLikeError >= 0.0)
      {
        //humanLike_legal = 0;
        return AutograspError;
      }else{
        //humanLike_legal = 1;
        std::cout<<"Human-like error1: "<<HumanLikeError<<std::endl;
        return HumanLikeError;
      }
    }else
    {
      if(HumanLikeError >= 0.0)
      {
        //humanLike_legal = 0;
        return ContactError;
      }else{
        //humanLike_legal = 1;
        std::cout<<"Human-like error2: "<<HumanLikeError<<std::endl;
        return ContactError;
      }
    }
  }else{
    return ContactError;
  }
}

/*! This version moves the palm in the direction of the object, attempting to establish contact on the palm
    before closing the fingers and establishing contacts on the finger.
*/
double
EBMGuidedAutoGraspQualityEnergy::approachAutograspQualityEnergy() const
{
  transf initialTran = mHand->getTran();
  bool contact = mHand->approachToContact(30);
  double autoQuality = 0.0;
  if (contact)
  {
    if (mHand->getPalm()->getNumContacts() == 0)
    {
      contact = false;
    }
  }
  if (!contact)
  {
    //if moving the hand in does not result in a palm contact, move out and grasp from the initial position
    //this allows us to obtain fingertip grasps if we want those
    DBGP("Approach found no contacts");
    mHand->setTran(initialTran);
  }
  else
  {
    DBGP("Approach results in contact");
  }
  autoQuality = autograspQualityEnergy();
  std::cout<<"Auto grasp energy: "<<autoQuality<<std::endl;

  return autoQuality;
}

/*! This function simply closes the hand and computes the real grasp quality that results
*/
double
EBMGuidedAutoGraspQualityEnergy::autograspQualityEnergy() const
{
  DBGP("Autograsp quality computation");
  mHand->autoGrasp(false, 1.0);
  mHand->getGrasp()->collectContacts();
  mHand->getGrasp()->updateWrenchSpaces();
  double volQual = mVolQual->evaluate();
  double epsQual = mEpsQual->evaluate();
  if (epsQual < 0)
  {
    epsQual = 0;
  } //returns -1 for non-FC grasps
  std::cout << "Autograsp quality: " << volQual << " volume and " << epsQual << " epsilon." << std::endl;
  return -(30 * volQual) - (100 * epsQual);
}

double
EBMGuidedAutoGraspQualityEnergy::ebm_pythonInterface(std::vector<double> &dofVals, int numDOF, std::string modelPath) const
{
  PyRun_SimpleString("import sys");
  //The file path of the ebmPythonInterface.py
  PyRun_SimpleString("sys.path.append('/root/WorkSpace/EBM_Hand')");
  pModule = PyImport_ImportModule("ebmPythonInterface");

  if (!pModule)
  {
    std::cout << "pModle is null" << std::endl;
  }

  if (PyImport_ImportModule("ebmPythonInterface") == NULL || PyErr_Occurred())
  {
    PyErr_Print();
  }

  while (!pModule)
  {
    sleep(3.0);
  }

  if (!pModule)
  {
    std::cout << "Error: python module is null!" << std::endl;
    return 0;
  }
  else
  {
    std::cout << "python module is successful" << std::endl;
  }

  pFunc = PyObject_GetAttrString(pModule, "dof_ebm");
  if (!pFunc)
  {
    std::cout << "Error: python pFunc_dof_ebm is null!" << std::endl;
    return 0;
  }
  else
  {
    std::cout << "python pFunc_dof_ebm is successful" << std::endl;
  }

  //创建参数:
  pArgs = PyTuple_New(2); //函数调用的参数传递均是以元组的形式打包的,2表示参数个数
  list_dof = PyList_New(0);

  for (int i = 0; i < numDOF; i++)
  {
    double degVal = dofVals[i] * 57.3;
    //std::cout<<degVal<<std::endl;
    PyList_Append(list_dof, Py_BuildValue("d", degVal));
  }

  PyTuple_SetItem(pArgs, 0, list_dof);

  PyTuple_SetItem(pArgs, 1, Py_BuildValue("s", modelPath.c_str()));

  //返回值
  pReturn = PyEval_CallObject(pFunc, pArgs); //调用函数
  //将返回值转换为double类型
  double result = 0.0;
  PyArg_Parse(pReturn, "d", &result); //d表示转换成double型变量
  std::cout << "ebm value: " << result << std::endl;

  Py_DECREF(pArgs);
  Py_DECREF(pModule);
  Py_DECREF(pFunc);
  Py_DECREF(list_dof);
  Py_DECREF(pReturn);

  pModule = NULL;
  pFunc = NULL;
  pReturn = NULL;
  pArgs = NULL;
  list_dof = NULL;
  return result;
}

double
EBMGuidedAutoGraspQualityEnergy::contactEnergy(int &closeContacts) const
{
  //DBGP("Contact energy computation")
  //average error per contact
  VirtualContact *contact;
  vec3 p, n, cn;
  double totalError = 0;
  mHand->getGrasp()->collectVirtualContacts();
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

    if (fabs(dist) < 20 && d < 0.3)
    {
      closeContacts++;
    }
  }

  totalError /= mHand->getGrasp()->getNumContacts();

  //DBGP("Contact energy: " << totalError);
  std::cout<<"Contact energy: "<<totalError<<std::endl;
  std::cout<<"Number of closed contacts: "<<closeContacts<<std::endl;
  return totalError;
}

double
EBMGuidedAutoGraspQualityEnergy::potentialQualityEnergy() const
{
  bool verbose = false;
  VirtualContact *contact;
  vec3 p, n, cn;
  int count = 0;
  //DBGP("Potential quality energy computation")
  for (int i = 0; i < mHand->getGrasp()->getNumContacts(); i++)
  {
    contact = (VirtualContact *)mHand->getGrasp()->getContact(i);
    contact->computeWrenches(true, false);
    contact->getObjectDistanceAndNormal(mObject, &p, NULL);
    n = contact->getWorldNormal();
    double dist = p.norm();
    p = p.normalized();
    double cosTheta = n.dot(p);
    double factor = potentialQualityScalingFunction(dist, cosTheta);
    if (verbose)
    {
      fprintf(stderr, "VC %d on finger %d link %d\n", i, contact->getFingerNum(), contact->getLinkNum());
      fprintf(stderr, "Distance %f cosTheta %f\n", dist, cosTheta);
      fprintf(stderr, "Scaling factor %f\n\n", factor);
    }
    contact->scaleWrenches(factor);
    if (factor > 0.25)
    {
      count++;
      contact->mark(true);
    }
    else
    {
      contact->mark(false);
    }
  }
  double gq = -1;
  //to make computations more efficient, we only use a 3D approximation
  //of the 6D wrench space
  std::vector<int> forceDimensions(6, 0);
  forceDimensions[0] = forceDimensions[1] = forceDimensions[2] = 1;
  if (count >= 3)
  {
    mHand->getGrasp()->updateWrenchSpaces(forceDimensions);
    gq = mEpsQual->evaluate();
  }
  if (verbose)
  {
    fprintf(stderr, "Quality: %f\n\n", gq);
  }
  if (count)
  {
    DBGP("Count: " << count << "; Gq: " << gq << ";");
  }
  return -gq;
  return 10;
}

double
EBMGuidedAutoGraspQualityEnergy::potentialQualityScalingFunction(double dist, double cosTheta) const
{
  double sf = 0;
  /*
  (void*)&cosTheta;
  if (dist<0) return 0;
  sf += 100.0 / ( pow(2.0,dist/15.0) );
  //comment this out if you don't care about normals
  //sf += 100 - (1 - cosTheta) * 100.0 / 2.0;
  */
  /*
  if (cosTheta < 0.8) return 0; //about 35 degrees
  sf = 10.0 / ( pow((double)3.0,dist/25.0) );
  if (sf < 0.25) sf = 0; //cut down on computation for tiny values. more than 50mm away
  return sf;
  */
  /*
  if (cosTheta < 0.8) return 0; //about 35 degrees
  if (dist < 20) sf = 1;
  else {
  sf = 1.5 - 1.5 / ( pow((double)3.0,(dist-20)/25.0) );
  sf = cos( sf*sf ) + 1;
  }
  sf = sf*10;
  if (sf < 0.25) sf = 0; //cut down on computation for tiny values. more than 50mm away
  return sf;
  */
  if (cosTheta < 0.7)
  {
    return 0;
  }
  if (dist > 50)
  {
    return 0;
  }
  sf = cos(3.14 * dist / 50.0) + 1;
  return sf;
}

double
EBMGuidedAutoGraspQualityEnergy::ebmEnergy() const
{
  //Evaluating human-like energy of the grasp candiate after closing the human hand
  double ebmQuality = 0.0;
  int numDof = mHand->getNumDOF();
  double *dofVals = new double[numDof];
  mHand->getDOFVals(dofVals);

  std::vector<double> DOFs(numDof);
  for (int i = 0; i < numDof; i++)
  {
    DOFs[i] = dofVals[i];
  }
  delete[] dofVals;
  ebmQuality = this->ebm_pythonInterface(DOFs, numDof, mHand->getEBMPath());
  return ebmQuality;
}