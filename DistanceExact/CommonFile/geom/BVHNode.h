#ifndef BVH_NODE_H
#define BVH_NODE_H

#include "../MathBasic.h"
#include "../IOBasic.h"

PRJ_BEGIN

template <typename T,typename BBOX=BBox<scalar> >
struct ALIGN_16 Node : public SerializableBase {
  static const int dim=BBOX::dim;
  using SerializableBase::read;
  using SerializableBase::write;
  using SerializableBase::operator<;
  using SerializableBase::operator=;
  typedef BBOX BoxType;
  EIGEN_DEVICE_FUNC Node():_l(-1),_r(-1),_parent(-1),_nrCell(-1) {}
  bool read(std::istream& is,IOData* dat) {
    readBinaryData(_bb,is);
    readBinaryData(_cell,is,dat);
    readBinaryData(_l,is);
    readBinaryData(_r,is);
    readBinaryData(_parent,is);
    readBinaryData(_nrCell,is);
    return is.good();
  }
  bool write(std::ostream& os,IOData* dat) const {
    writeBinaryData(_bb,os);
    writeBinaryData(_cell,os,dat);
    writeBinaryData(_l,os);
    writeBinaryData(_r,os);
    writeBinaryData(_parent,os);
    writeBinaryData(_nrCell,os);
    return os.good();
  }
  std::shared_ptr<SerializableBase> copy() const {
    return std::shared_ptr<SerializableBase>(new Node<T,BBOX>);
  }
  virtual std::string type() const override {
    return typeid(Node<T,BBOX>).name();
  }
  EIGEN_DEVICE_FUNC Node<T,BBOX>& operator=(const Node<T,BBOX>& other) {
    _bb=other._bb;
    _cell=other._cell;
    _l=other._l;
    _r=other._r;
    _parent=other._parent;
    _nrCell=other._nrCell;
    return *this;
  }
  ALIGN_16 BBOX _bb;
  ALIGN_16 T _cell;
  ALIGN_16 sizeType _l,_r,_parent,_nrCell;
};

PRJ_END

#endif
