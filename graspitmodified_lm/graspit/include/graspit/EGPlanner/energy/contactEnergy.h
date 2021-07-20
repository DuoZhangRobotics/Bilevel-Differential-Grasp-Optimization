
#ifndef _contactenergy_h_
#define _contactenergy_h_

#include "graspit/EGPlanner/energy/searchEnergy.h"
#include <string.h>
#include <vector>

class ContactEnergy: public SearchEnergy
{
  public:
    double energy() const;
  /**
   * @brief An interface to run EBM model.
   * 
   * @param dofVals DOF of hand
   * @param numDOF number of DOF
   * @param modelPath The file path of ebm model
   * @return double ebm energy value that is used to predict human-like grasping.
   */
  //  double method1(std::vector<double>& dofVals, int numDOF, std::string modelPath) const;
};


#endif
