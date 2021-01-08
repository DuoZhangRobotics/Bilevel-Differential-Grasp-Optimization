#include <TrajOpt/PhaseSequence.h>
//Position
#include <TrajOpt/PositionSequence.h>
#include <TrajOpt/SampledPositionSequence.h>
#include <TrajOpt/ConstantPositionSequence.h>
#include <TrajOpt/PenetratingPositionSequence.h>
#include <TrajOpt/PositionSequenceAboveGround.h>
#include <TrajOpt/PositionSequenceDense.h>
//Force
#include <TrajOpt/ForceSequence.h>
#include <TrajOpt/PBTOForceSequence.h>
#include <TrajOpt/PBTOForceSequencePose.h>
#include <TrajOpt/ConstantForceSequence.h>
#include <TrajOpt/SampledForceSequence.h>
#include <TrajOpt/ForceSequenceInCone.h>
#include <TrajOpt/ForceSequenceKKT.h>
#include <TrajOpt/LPForceSequence.h>
#include <TrajOpt/LPForceSequenceMollifier.h>
#include <TrajOpt/QPForceSequenceMollifier.h>
//Objective
#include <TrajOpt/PositionSequenceDebugObj.h>
#include <TrajOpt/ForceSequenceDebugObj.h>
#include <TrajOpt/PositionObj.h>
#include <TrajOpt/VelocityObj.h>
#include <TrajOpt/ShuffleAvoidanceObj.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

#define FORCE_CREATE true
template <typename T,int DIM>
void debug(sizeType trial,T thres)
{
  std::shared_ptr<EnvironmentHeight<T>> env;
  if(exists("envHeight.dat") && !FORCE_CREATE) {
    env.reset(new EnvironmentHeight<T>(0.1f));
    env->SerializableBase::read("envHeight.dat");
  } else {
    env.reset(new EnvironmentHeight<T>(0.1f));
    //env->createStair(5,5,5,2.5,0.1,10);
    env->createHills(5,5,[&](scalarD x,scalarD y) {
      return std::sin(x)*std::sin(y);
    },16);
    env->SerializableBase::write("envHeight.dat");
  }

  T mu=0.75f;
  T dt=0.05f;
  T phi0=0.01f;
  typename ForceSequenceDebugObj<T,DIM>::Vec3T dir;
  for(sizeType fixedSpan=0; fixedSpan<2; fixedSpan++)
    for(sizeType even=0; even<2; even++) {
      DSSQPObjectiveCompound<T> obj;
      std::shared_ptr<PhaseSequence<T>> seq(new PhaseSequence<T>(obj,"",0.5,5.0,fixedSpan?7:-7,4));
      if(!seq->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      //position
      std::shared_ptr<PositionSequence<T,DIM>> posSeq(new PositionSequence<T,DIM>(obj,"",phi0,even,seq,env));
      if(!posSeq->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      std::shared_ptr<ConstantPositionSequence<T,DIM>> cposSeq(new ConstantPositionSequence<T,DIM>(obj,"",phi0,true,seq,env));
      if(!cposSeq->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      std::shared_ptr<SampledPositionSequence<T,DIM>> sposSeq(new SampledPositionSequence<T,DIM>(obj,"",phi0,dt,seq,env));
      if(!sposSeq->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      std::shared_ptr<PenetratingPositionSequence<T,DIM>> pposSeq(new PenetratingPositionSequence<T,DIM>(obj,"",phi0,seq,env));
      if(!pposSeq->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      std::shared_ptr<PositionSequenceAboveGround<T,DIM>> posSeqAboveGround(new PositionSequenceAboveGround<T,DIM>(obj,"",true,dt,posSeq));
      if(!posSeqAboveGround->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      std::shared_ptr<PositionSequenceDense<T,DIM>> posSeqDense(new PositionSequenceDense<T,DIM>(obj,"",dt,phi0,posSeq));
      if(!posSeqDense->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      dir.setRandom();
      std::shared_ptr<PositionObj<T,DIM>> posObj(new PositionObj<T,DIM>(obj,"",dt,dir,1,posSeq));
      if(!posObj->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      posObj.reset(new PositionObj<T,DIM>(obj,"",dt,dir,0,posSeq));
      if(!posObj->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      dir.setRandom();
      std::shared_ptr<VelocityObj<T,DIM>> velObj(new VelocityObj<T,DIM>(obj,"",dt*2,dt,dir,1,posSeq));
      if(!velObj->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      velObj.reset(new VelocityObj<T,DIM>(obj,"",dt*2,dt,dir,0,posSeq));
      if(!velObj->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      //force
      std::shared_ptr<ForceSequence<T,DIM>> forceSeq(new ForceSequence<T,DIM>(obj,"",mu,!even,seq,env));
      if(!forceSeq->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      std::shared_ptr<PBTOForceSequence<T,DIM>> forceSeqPBTO(new PBTOForceSequence<T,DIM>(obj,"",dt,100,100,posSeq));
      if(!forceSeqPBTO->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      {
        sizeType nrJ=10;
        typename C2EnvWrenchConstructor<T>::Vec3T G(0,0,-9.81f);
        std::shared_ptr<ArticulatedBody> body(new ArticulatedBody());
        body->randomize(nrJ);
        std::vector<std::shared_ptr<PositionSequence<T,DIM>>> poses;
        poses.push_back(posSeq);
        std::vector<std::shared_ptr<ForceSequence<T,DIM>>> forces;
        forces.push_back(forceSeq);
        std::shared_ptr<C2EnvWrenchConstructor<T>> wrench(new C2EnvWrenchConstructor<T>(*body,env,G));
        wrench->_externalForces.push_back(EndEffectorBounds(nrJ-1,Vec3d::Random(),0.1f));
        std::shared_ptr<PBDDynamicsSequence<T,DIM>> PBD(new PBDDynamicsSequence<T,DIM>(obj,"",dt,DSSQPObjectiveComponent<T>::infty(),body,G,poses,forces,wrench->_externalForces));
        {
          std::shared_ptr<LPForceSequence<T,DIM>> forceSeqLP(new LPForceSequence<T,DIM>(obj,"",wrench,posSeq));
          forceSeqLP->setPBDDynamicsSequence(PBD);
          if(!forceSeqLP->debug(obj.inputs(),trial,thres)) {
            ASSERT(false)
          }
        }
        {
          std::shared_ptr<LPForceSequenceMollifier<T,DIM>> forceSeqLP(new LPForceSequenceMollifier<T,DIM>(obj,"",wrench,forceSeq));
          forceSeqLP->setPBDDynamicsSequence(PBD);
          if(!forceSeqLP->debug(obj.inputs(),trial,thres)) {
            ASSERT(false)
          }
        }
        {
          std::shared_ptr<QPForceSequenceMollifier<T,DIM>> forceSeqQP(new QPForceSequenceMollifier<T,DIM>(obj,"",wrench,forceSeq));
          forceSeqQP->setPBDDynamicsSequence(PBD);
          if(!forceSeqQP->debug(obj.inputs(),trial,thres)) {
            ASSERT(false)
          }
        }
        std::shared_ptr<PBTOForceSequencePose<T,DIM>> forceSeqPBTOPose(new PBTOForceSequencePose<T,DIM>(obj,"",dt,100,100,posSeq,wrench->_externalForces.back()));
        forceSeqPBTOPose->setPBDDynamicsSequence(PBD);
        if(!forceSeqPBTOPose->debug(obj.inputs(),trial,thres)) {
          ASSERT(false)
        }
      }
      std::shared_ptr<ConstantForceSequence<T,DIM>> cforceSeq(new ConstantForceSequence<T,DIM>(obj,"",mu,seq,env));
      if(!cforceSeq->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      std::shared_ptr<SampledForceSequence<T,DIM>> sforceSeq(new SampledForceSequence<T,DIM>(obj,"",mu,dt,6,SampledForceSequence<T,DIM>::Vec3T::UnitZ(),posSeq,env));
      if(!sforceSeq->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      std::shared_ptr<ForceSequenceInCone<T,DIM>> forSeqInCone(new ForceSequenceInCone<T,DIM>(obj,"",mu,posSeq,forceSeq));
      if(!forSeqInCone->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      std::shared_ptr<ForceSequenceKKT<T,DIM>> forSeqKKT(new ForceSequenceKKT<T,DIM>(obj,"",sposSeq,sforceSeq));
      if(!forSeqKKT->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      //objective
      std::shared_ptr<PositionSequenceDebugObj<T,DIM>> posSeqDebugObj(new PositionSequenceDebugObj<T,DIM>(obj,"",dt,1,posSeq));
      if(!posSeqDebugObj->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      std::shared_ptr<ShuffleAvoidanceObj<T,DIM>> shuffleAoidanceObj(new ShuffleAvoidanceObj<T,DIM>(obj,"",dt,1,posSeq));
      if(!shuffleAoidanceObj->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
      dir.setRandom();
      dir/=std::sqrt(dir.squaredNorm());
      std::shared_ptr<ForceSequenceDebugObj<T,DIM>> forceSeqDebugObj(new ForceSequenceDebugObj<T,DIM>(obj,"",dt,dir*M_PI*2,1,forceSeq));
      if(!forceSeqDebugObj->debug(obj.inputs(),trial,thres)) {
        ASSERT(false)
      }
    }
}
int main()
{
  OmpSettings::getOmpSettingsNonConst().setNrThreads(1);
  mpfr_set_default_prec(1024U);
  debug<mpfr::mpreal,2>(10,0     );
  debug<mpfr::mpreal,3>(10,0     );
  //debug<mpfr::mpreal,2>(10,-1E-4f);
  //debug<mpfr::mpreal,3>(10,-1E-4f);
  return 0;
}
