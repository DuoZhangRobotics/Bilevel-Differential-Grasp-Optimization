//######################################################################
//
// GraspIt!
// Copyright (C) 2002-2009  Columbia University in the City of New York.
// All rights reserved.
//
// GraspIt! is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GraspIt! is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GraspIt!.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s): Andrew T. Miller
//
// $Id: main.cpp,v 1.16 2010/01/13 23:09:30 cmatei Exp $
//
//######################################################################

/*! \mainpage GraspIt! Developer Documentation
  \image html logo.jpg

  These pages document the GraspIt! source code. Please remember this is
  research code. There are still plenty of pieces of code that are unfinished
  and several bugs that need to be fixed.

  More information and original source code for the packages included with
  GraspIt! can be found in the following places:

  - <b>qhull:</b> http://www.qhull.org
  - <b>maxdet:</b> http://www.stanford.edu/~boyd/old_software/MAXDET.html
  - <b>tinyxml:</b> http://www.grinninglizard.com/tinyxml/
*/

/*! \file
  \brief Program execution starts here.  Server is started, main window is built,
  and the interactive loop is started.

  The main call returns an exit code as indicated by the graspit GUI (0 by default)
  to provide feedback to calling program, if desired.
 */

#include <iostream>
#include "graspit/world.h"
#include "graspit/graspitCore.h"
#include "graspit/graspitParser.h"
#include "graspit/body.h"
#include "graspit/robot.h"
#include "graspit/EGPlanner/egPlanner.h"
#include "graspit/EGPlanner/searchState.h"
#include "graspit/EGPlanner/simAnnPlanner.h"
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <strstream>
#include "graspit/quality/qualEpsilon.h"
#include "graspit/quality/quality.h"
#include "graspit/grasp.h"
#include <fstream>
#include <string>
class SimAnnPlannerCMD : public SimAnnPlanner
{
public:
  SimAnnPlannerCMD(Hand* h)
    :SimAnnPlanner(h) {}
  bool runCMD() {
    bool done=false;
    while(!done) {
      PlannerState s=getState();
      switch(s) {
      case STARTING_THREAD: //do nothing
        break;
      case INIT:
        sleep(0.1);
        break;
      case READY:
        sleep(0.1);
        break;
      case RUNNING:
        mainLoop();
        break;
      case DONE:
        done = true;
        break;
      case EXITED: //Do nothing
        break;
      }
      if(!done && checkTerminationConditions())
        break;
      if(getCurrentStep()%200 == 0)
        std::cout << "Run " << getCurrentStep() << "/" << mMaxSteps << std::endl;
    }
    setState(EXITED);
    return done;
  }
};
transf loadTransf(const cmdline::parser& parser,const std::string& rot,const std::string& trans)
{
  transf tf;
  Quaternion q=Quaternion::Identity();
  vec3 t=vec3::Zero();
  if(parser.exist(rot)) {
    std::string val=parser.get<std::string>(rot);
    std::replace(val.begin(),val.end(),',',' ');
    std::istringstream(val) >> q.w() >> q.x() >> q.y() >> q.z();
  }
  if(parser.exist(trans)) {
    std::string val=parser.get<std::string>(trans);
    std::replace(val.begin(),val.end(),',',' ');
    std::istringstream(val) >> t.x() >> t.y() >> t.z();
  }
  tf.set(q,t);
  return tf;
}
void loadBodyFile(GraspitCore& core,int argc,char **argv)
{
  cmdline::parser parsed_args;
  parsed_args.add<std::string>("bodyFile");
  parsed_args.add<std::string>("bodyRot");
  parsed_args.add<std::string>("bodyTrans");
  parsed_args.parse(argc,argv);
  World* world=core.getWorld();
  if(parsed_args.exist("bodyFile")) {
    std::string file=parsed_args.get<std::string>("bodyFile");
    if(file.length() > 4) {
      Body* body=world->importBody("GraspableBody",QString(file.c_str()));
      if(!body) {
        std::cout << "Cannot load bodyFile: " << file << std::endl;
        exit(EXIT_FAILURE);
      }
      transf bodyTrans=loadTransf(parsed_args,"bodyRot","bodyTrans");
      body->setTran(bodyTrans);
    } else {
      std::cout << "Invalid bodyFile!" << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    std::cout << "Cannot find bodyFile!" << std::endl;
    exit(EXIT_FAILURE);
  }
}
void loadRobotFile(GraspitCore& core,int argc,char **argv)
{
  cmdline::parser parsed_args;
  parsed_args.add<std::string>("robotFile");
  parsed_args.add<std::string>("robotRot");
  parsed_args.add<std::string>("robotTrans");
  parsed_args.add<std::string>("robotDOF");
  parsed_args.parse(argc,argv);
  World* world=core.getWorld();
  if(parsed_args.exist("robotFile")) {
    std::string file=parsed_args.get<std::string>("robotFile");
    if(file.length() > 4) {
      QString pathQ(file.c_str());
      Robot* robot=world->importRobot(pathQ);
      if(!robot) {
        std::cout << "Cannot load robotFile: " << file << std::endl;
        exit(EXIT_FAILURE);
      }
      transf robotTrans=loadTransf(parsed_args,"robotRot","robotTrans");
      robot->setTran(robotTrans);
      if(parsed_args.exist("robotDOF")) {
        std::string val=parsed_args.get<std::string>("robotDOF");
        std::replace(val.begin(),val.end(),',',' ');
        QString valQ(val.c_str());
        QTextStream is(&valQ);
        robot->readDOFVals(is);
      }
    } else {
      std::cout << "Invalid robotFile!" << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    std::cout << "Cannot find robotFile!" << std::endl;
    exit(EXIT_FAILURE);
  }
}


std::vector<string> &split(const string &str, char delim, std::vector<string> &elems, bool skip_empty = true) {
    std::istringstream iss(str);
    for (string item; std::getline(iss, item, delim); )
        if (skip_empty && item.empty()) continue;
        else elems.push_back(item);
    return elems;
}

double * desiredDOF(int argc,char **argv, int numDOF)
{
  double * desiredvals = new double[numDOF];
  cmdline::parser parsed_args;
  parsed_args.add<std::string>("desiredDOF");
  parsed_args.parse(argc,argv);
  if(parsed_args.exist("desiredDOF"))
  {
    std::string val=parsed_args.get<std::string>("desiredDOF");
    std::vector<string> result;
    split(val, ',', result);
    for (int i = 0; i < result.size(); ++i)
    {
        std::strstream ss;
        ss<<result[i];
        ss>>desiredvals[i];
        std::cout<<desiredvals[i]<<std::endl;
    }


//    std::cout<<'is'<<std::endl;
  }
  else
  {
    std::cout << "Cannot find desiredDOF!" << std::endl;
    exit(EXIT_FAILURE);
  }
  return desiredvals;
}



int countContacts(Hand* h)
{
  int nrCT=0;
  std::list<Contact*> cc=h->getContacts();
  for(std::list<Contact*>::const_iterator beg=cc.begin(),end=cc.end(); beg!=end; beg++) {
    Contact* ct=*beg;
    if(dynamic_cast<GraspableBody*>(ct->getBody1()) && !dynamic_cast<GraspableBody*>(ct->getBody2()))
      nrCT++;
    else if(dynamic_cast<GraspableBody*>(ct->getBody2()) && !dynamic_cast<GraspableBody*>(ct->getBody1()))
      nrCT++;
  }
  return nrCT;
}
typedef char* ARG;
std::string output_path;
int main(int argc, char **argv)
{
  cmdline::parser parsed_args;
  parsed_args.add<std::string>("resultFile");
  parsed_args.add<int>("maxStep",0,"",true,40000);
  parsed_args.add<std::string>("EGEnergy",0,"",true,"STRICT_AUTO_GRASP_ENERGY");
  parsed_args.parse(argc,argv);

  std::cout<<parsed_args.get<std::string>("resultFile").c_str();
 
   //std::string output_path;
  output_path.assign(parsed_args.get<std::string>("resultFile").c_str());

  std::cout<<output_path;


  //create core
  int argc_tmp=2;
  ARG argv_tmp[2]= {
    (char*)"graspit",
    (char*)"--headless"
  };
  GraspitCore core(argc_tmp,argv_tmp);
  //load
  loadBodyFile(core,argc,argv);
  loadRobotFile(core,argc,argv);

  int numDOF = 18;
  double *desiredvals = desiredDOF(argc,argv,numDOF);

  //planner
  if(core.getWorld()->getNumGB() == 0) {
    std::cout << "World has no graspable objects!" << std::endl;
    exit(EXIT_FAILURE);
  }
  //state
  GraspPlanningState state(core.getWorld()->getCurrentHand());
  state.setPositionType(SPACE_AXIS_ANGLE);
  state.setObject(core.getWorld()->getGB(0));
  state.setRefTran(core.getWorld()->getGB(0)->getTran());
  state.reset();
  //planner
  SimAnnPlannerCMD planner(core.getWorld()->getCurrentHand());
  planner.setModelState(&state);
  std::string energy=parsed_args.get<std::string>("EGEnergy");
  if(energy == "CONTACT_ENERGY" || energy == "POTENTIAL_QUALITY_ENERGY" ||
     energy == "AUTO_GRASP_QUALITY_ENERGY" || energy == "GUIDED_POTENTIAL_QUALITY_ENERGY" ||
     energy == "GUIDED_AUTO_GRASP_QUALITY_ENERGY" || energy == "STRICT_AUTO_GRASP_ENERGY" ||
     energy == "COMPLIANT_ENERGY" || energy == "DYNAMIC_AUTO_GRASP_ENERGY")
    planner.setEnergyType(energy);
  else {
    std::cout << "Unknown energy type: " << energy << std::endl;
    exit(EXIT_FAILURE);
  }
  planner.setContactType(CONTACT_LIVE);
  planner.setMaxSteps(parsed_args.get<int>("maxStep"));
  if(!planner.resetPlanner()) {
    std::cout << "Error reset EGPlanner!" << std::endl;
    exit(EXIT_FAILURE);
  }

  planner.startPlanner();
  core.getWorld()->getCurrentHand()->auto_target(false,1,false,desiredvals);
  core.getWorld()->getCurrentHand()->autoGrasp(false);

  Hand *h;
  Grasp *g;
  h=core.getWorld()->getHand(0);
  g=h->getGrasp();
  QualEpsilon mEpsQual(g, QString("LM"), "L1 Norm");
  g->update();
  double val;
  val=mEpsQual.evaluate();
  std::cout<<val<<std::endl;

  if(parsed_args.exist("resultFile")) {
    std::ostringstream oss;
    oss << parsed_args.get<std::string>("resultFile").c_str()<< "grasp" << ".xml";
    core.getWorld()->save(QString(oss.str().c_str()));
  }
  
  std::ofstream myfile;
  std::string myString(parsed_args.get<std::string>("resultFile").c_str());
  myString.append("eps_qual.txt");
  myfile.open(myString.c_str());
  std::ostringstream os;
  os << val;
  std::string str = os.str();
  myfile<<str.c_str();
  myfile.close();
  
  exit(EXIT_SUCCESS);
  return 0;
}
