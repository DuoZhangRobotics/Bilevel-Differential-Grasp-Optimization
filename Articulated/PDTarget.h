#ifndef PD_TARGET_H
#define PD_TARGET_H

#include <CommonFile/IO.h>
#include <Utils/SparseUtils.h>

PRJ_BEGIN

struct ArticulatedBody;
#include <Utils/ArticulatedBodyPragma.h>
struct PDTarget : public SerializableBase
{
public:
  typedef scalarD T;
  DECL_MAP_TYPES_T
  PDTarget();
  PDTarget(const Vec& PCoef,const Vec& DCoef,const std::vector<std::tuple<scalarD,Vec,Vec>>& spline);
  PDTarget(const Vec& PCoef,const Vec& DCoef,const Vec& s);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  void writeVTKSeq(const ArticulatedBody& body,const std::string& path,T dt);
  void reset();
  void advance(T dt);
  Vec s() const;
  void setS(const Vec& s);
  Vec splineInterp() const;
  Vec _PCoef,_DCoef;
private:
  void advanceSpline();
  std::vector<std::tuple<T,Vec,Vec>> _spline;
  sizeType _id;
  Vec _s;
  T _t;
};

PRJ_END

#endif
