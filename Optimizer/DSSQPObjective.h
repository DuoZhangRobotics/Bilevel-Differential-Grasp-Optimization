#ifndef DSSQP_OBJECTIVE_H
#define DSSQP_OBJECTIVE_H

#include <Utils/ArticulatedBodyPragma.h>
#include <Utils/SparseUtils.h>
#include <unordered_map>

PRJ_BEGIN

template <typename T>
struct DSSQPObjective;
template <typename T>
struct DSSQPVariable;
template <typename T>
struct DSSQPObjectiveComponent;
template <typename T>
struct DSSQPObjectiveCompound;

//definition
template <typename T>
struct DSSQPObjective : public SparseTraits<T>
{
public:
  using typename SparseTraits<T>::Vec;
  using typename SparseTraits<T>::DMat;
  using typename SparseTraits<T>::SMat;
  using typename SparseTraits<T>::STrip;
  using typename SparseTraits<T>::STrips;
  static T infty();
  virtual ~DSSQPObjective() {}
  virtual Vec lb() const;
  virtual Vec ub() const;
  virtual Vec gl() const;
  virtual Vec gu() const;
  virtual Vec init() const;
  int nnzJ();
  //constraint
  virtual int operator()(const Vec& x,Vec& fvec,DMat* fjac);
  virtual int operator()(const Vec& x,Vec& fvec,SMat* fjac);
  virtual int operator()(const Vec& x,Vec& fvec,STrips* fjac);
  //objective
  virtual T operator()(const Vec& x,Vec* fgrad);
  //problem size
  virtual int inputs() const;
  virtual int values() const;
protected:
  SMat _tmpFjacs;
  STrips _tmp;
};
template <typename T>
struct DSSQPVariable
{
  sizeType _id;
  T _l,_u,_init;
};
template <typename T>
struct DSSQPObjectiveComponent : public DSSQPObjective<T>
{
public:
  DECL_MAP_TYPES_T
  typedef std::unordered_map<sizeType,std::string> VARMAPINV;
  typedef std::unordered_map<std::string,DSSQPVariable<T>> VARMAP;
  typedef std::unordered_map<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>> CONSMAP;
  DSSQPObjectiveComponent(DSSQPObjectiveCompound<T>& obj,const std::string& name,bool force=false);
  virtual ~DSSQPObjectiveComponent();
  virtual int inputs() const override;
  virtual int values() const override;
  //constraint
  virtual int operator()(const Vec& x,Vec& fvec,STrips* fjac) override;
  //objective
  virtual T operator()(const Vec& x,Vec* fgrad) override;
  virtual bool debug(sizeType inputs,sizeType nrTrial,T thres);
  virtual Vec makeValid(const Vec& x) const;
  virtual bool debugGradientConditional(const std::string& entryStr,T ref,T err,T thres) const;
  virtual bool debugGradientConditionalVID(const std::string& entryStr,sizeType vid,T ref,T err,T thres) const;
  void setOffset(sizeType offset);
  std::string _name;
  std::vector<T> _gl,_gu;
  const VARMAP* _vars;
  sizeType _offset;
private:
  sizeType _inputs;
};
template <typename T>
struct DSSQPObjectiveCompound : public DSSQPObjective<T>
{
  enum VARIABLE_OP
  {
    MUST_NEW,
    MUST_EXIST,
    NEW_OR_EXIST,
  };
  using typename SparseTraits<T>::Vec;
  using typename SparseTraits<T>::DMat;
  using typename SparseTraits<T>::SMat;
  using typename SparseTraits<T>::STrip;
  using typename SparseTraits<T>::STrips;
  typedef std::unordered_map<sizeType,std::string> VARMAPINV;
  typedef std::unordered_map<std::string,DSSQPVariable<T>> VARMAP;
  typedef std::unordered_map<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>> CONSMAP;
  virtual ~DSSQPObjectiveCompound() {}
  virtual Vec lb() const override;
  virtual Vec ub() const override;
  virtual Vec gl() const override;
  virtual Vec gu() const override;
  virtual Vec init() const override;
  virtual int inputs() const override;
  virtual int values() const override;
  //constraint
  virtual int operator()(const Vec& x,Vec& fvec,STrips* fjac) override;
  //objective
  virtual T operator()(const Vec& x,Vec* fgrad) override;
  const DSSQPVariable<T>& addVar(const std::string& name,T l,T u,VARIABLE_OP op=NEW_OR_EXIST);
  const DSSQPVariable<T>& addVar(const std::string& name,VARIABLE_OP op=NEW_OR_EXIST);
  const DSSQPVariable<T>& addVar(sizeType id) const;
  DSSQPVariable<T>& addVar(sizeType id);
  void setVarInit(const std::string& name,T init);
  void setVarInit(sizeType id,T init);
  void checkViolation(const Vec* at=NULL);
  const CONSMAP& components() const;
  void addComponent(std::shared_ptr<DSSQPObjectiveComponent<T>> c);
  virtual bool debug(const std::string& str,sizeType nrTrial,T thres);
protected:
  CONSMAP _components;
  VARMAPINV _varsInv;
  VARMAP _vars;
};

PRJ_END

#endif
