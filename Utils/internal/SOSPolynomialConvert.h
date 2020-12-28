//fromString
template <typename T>
void convertFormatted(const std::string& str,T& val) {
  val=std::stod(str);
}
template <typename T>
void convertFormatted(const std::string& str,COMMON::SMatNumber<T>& val) {
  ASSERT_MSG(false,"Converting SMatNumber from string is not supported!")
}
template <typename T,char LABEL>
void convertFormatted(const std::string& str,SOSTerm<T,LABEL>& val) {
  val.readFormattedString(str);
}
template <typename T,char LABEL>
void convertFormatted(const std::string& str,SOSPolynomial<T,LABEL>& val) {
  val.readFormattedString(str);
}
//fromString
template <typename T>
void convert(const std::string& str,T& val) {
  val=std::stod(str);
}
template <typename T>
void convert(const std::string& str,COMMON::SMatNumber<T>& val) {
  ASSERT_MSG(false,"Converting SMatNumber from string is not supported!")
}
template <typename T,char LABEL>
void convert(const std::string& str,SOSTerm<T,LABEL>& val) {
  val << str;
}
template <typename T,char LABEL>
void convert(const std::string& str,SOSPolynomial<T,LABEL>& val) {
  val << str;
}

//toString
template <typename T>
std::string convertFormatted(T val) {
  return std::to_string(val);
}
template <typename T>
std::string convertFormatted(const COMMON::SMatNumber<T>& val) {
  ASSERT_MSG(false,"Converting SMatNumber to string is not supported!")
  return std::string();
}
template <typename T,char LABEL>
std::string convertFormatted(const SOSTerm<T,LABEL>& val) {
  return val.formattedString();
}
template <typename T,char LABEL>
std::string convertFormatted(const SOSPolynomial<T,LABEL>& val) {
  return val.formattedString();
}
//toString
template <typename T>
std::string convert(T val,bool& bracketNeeded) {
  bracketNeeded=false;
  return std::to_string(val);
}
template <typename T>
std::string convert(const COMMON::SMatNumber<T>& val) {
  ASSERT_MSG(false,"Converting SMatNumber to string is not supported!")
  return std::string();
}
template <typename T,char LABEL>
std::string convert(const SOSTerm<T,LABEL>& val,bool& bracketNeeded) {
  bracketNeeded=false;
  return val.toString();
}
template <typename T,char LABEL>
std::string convert(const SOSPolynomial<T,LABEL>& val,bool& bracketNeeded) {
  bracketNeeded=(sizeType)val._terms.size()>1;
  return val.toString();
}
