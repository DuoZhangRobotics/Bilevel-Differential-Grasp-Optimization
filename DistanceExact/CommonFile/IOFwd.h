#ifndef IO_FWD_H
#define IO_FWD_H

#include "MathBasic.h"
#include <memory>

PRJ_BEGIN

struct IOData;
struct ALIGN_16 SerializableBase {
  EIGEN_DEVICE_FUNC explicit SerializableBase() {}
  EIGEN_DEVICE_FUNC virtual ~SerializableBase() {}
  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  virtual bool read(std::istream& is,IOData* dat);
  virtual bool write(std::ostream& os,IOData* dat) const;
  virtual bool read(const std::string& path);
  virtual bool write(const std::string& path) const;
  virtual std::shared_ptr<SerializableBase> copy() const;
  virtual std::string type() const=0;
  EIGEN_DEVICE_FUNC virtual bool operator<(const SerializableBase& other) const {
    return type()<other.type();
  }
};
struct ALIGN_16 Serializable : public SerializableBase {
  using SerializableBase::read;
  using SerializableBase::write;
  using SerializableBase::copy;
  using SerializableBase::type;
  EIGEN_DEVICE_FUNC explicit Serializable(const std::string& className):_type(className) {}
  virtual std::string type() const override;
  EIGEN_DEVICE_FUNC virtual Serializable& operator=(const Serializable& other) {
    return *this;
  }
protected:
  void setType(const std::string& type);
  void setSerializableType(const std::string& type);
private:
  ALIGN_16 std::string _type;
};
template<typename T>
struct VTKWriter;

PRJ_END

#endif
