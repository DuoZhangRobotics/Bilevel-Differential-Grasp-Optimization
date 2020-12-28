#include <stack>
#include <string>
#include <functional>

template <typename T,char LABEL>
struct SOSPolynomialLexer
{
  enum TOKEN
  {
    LB,
    RB,
    MUL,
    POW,
    ADD,
    SUB,
    FLT,
    INT,
    VAR,
  };
  typedef SOSPolynomial<T,LABEL> POLY;
  struct TestLexer
  {
    bool operator()(TOKEN t,const std::string& str) {
      std::cout << " " << str;
      return true;
    }
  };
  POLY parse(const std::string& str) {
    while(!_opStack.empty())
      _opStack.pop();
    while(!_pStack.empty())
      _pStack.pop();
    bool succ=lexer(str,[&](TOKEN t,const std::string& str) {
      return parseOp(t,str);
    });
    ASSERT_MSG(succ,"SOSPolynomial parser failed!")
    while(!_opStack.empty())
      if(!applyOp())
        return POLY();
    return _pStack.top();
  }
  void testLex(const std::string& str) {
    lexer(str,TestLexer());
    std::cout << std::endl;
  }
  static POLY pow(const POLY& base,const POLY& exp) {
    ASSERT(base.orderAll()>0)  //we do not allow you to have a pow of number (this is not a calculator)
    POLY ret=base;
    sizeType expInt=ScalarOfT<POLY>::template convertBack<sizeType>(exp);
    ASSERT(expInt>0)
    for(; expInt>1; expInt--)
      ret*=base;
    return ret;
  }
  bool lexer(const std::string& str,std::function<bool(TOKEN,const std::string&)> f)
  {
    bool neg=false;
    for(size_t i=0; i<str.size(); i++) {
      size_t last=i;
      if(str[i]=='(') {
        if(!f(LB,str.substr(i,1)))
          return false;
      }  else if(str[i]==')') {
        if(!f(RB,str.substr(i,1)))
          return false;
      } else if(str[i]=='*') {
        if(!f(MUL,str.substr(i,1)))
          return false;
      } else if(str[i]=='^') {
        if(!f(POW,str.substr(i,1)))
          return false;
      } else if(str[i]=='+') {
        if(!f(ADD,str.substr(i,1)))
          return false;
      } else if(str[i]=='-') {
        if(i==0 || str[i-1]=='(' || str[i-1]=='*' || str[i-1]=='^' || str[i-1]=='+')
          neg=true;
        else if(!f(SUB,str.substr(i,1)))
          return false;
      } else if((str[i]>='A'&&str[i]<='Z') || (str[i]>='a'&&str[i]<='z')) {
        if(neg) {
          last--;
          neg=false;
        }
        while(i+1<str.size() && str[i+1]>='0' && str[i+1]<='9')
          i++;
        if(!f(VAR,str.substr(last,i-last+1)))
          return false;
      } else if((str[i]>='0' && str[i]<='9')||str[i]=='.') {
        if(neg) {
          last--;
          neg=false;
        }
        bool dot=str[i]=='.';
        while(i+1<str.size() && ((str[i+1]>='0' && str[i+1]<='9')||str[i+1]=='.')) {
          if(str[i+1]=='.') {
            if(dot)
              return false;
            else dot=true;
          }
          i++;
        }
        if(!f(dot?FLT:INT,str.substr(last,i-last+1)))
          return false;
      } else return false;
    }
    return true;
  }
  bool parseOp(TOKEN t,const std::string& str) {
    if(t==MUL) {
      while(!_opStack.empty() && _opStack.top()==POW)
        if(!applyOp())
          return false;
      _opStack.push(t);
    } else if(t==POW) {
      _opStack.push(t);
    } else if(t==ADD||t==SUB) {
      while(!_opStack.empty() && (_opStack.top()==POW||_opStack.top()==MUL||_opStack.top()==SUB))
        if(!applyOp())
          return false;
      _opStack.push(t);
    } else if(t==LB)
      _opStack.push(t);
    else if(t==RB) {
      while(_opStack.top()!=LB)
        if(!applyOp())
          return false;
      _opStack.pop();
    } else if(t==FLT)
      _pStack.push(ScalarOfT<POLY>::convert(atof(str.c_str())));
    else if(t==INT)
      _pStack.push(ScalarOfT<POLY>::convert(atoi(str.c_str())));
    else if(t==VAR) {
      if(str[0]=='-')   //Do not say -a1. Instead, say -1*a1
        return false;
      _pStack.push(ScalarOfT<POLY>::convertVar(atoi(str.substr(1).c_str()),str[0]));
    } else return false;
    return true;
  }
  bool applyOp()
  {
    if(_pStack.size()<2||_opStack.size()<1)
      return false;
    POLY R=_pStack.top();
    _pStack.pop();
    POLY L=_pStack.top();
    _pStack.pop();
    if(_opStack.top()==ADD)
      _pStack.push(L+R);
    else if(_opStack.top()==SUB)
      _pStack.push(L-R);
    else if(_opStack.top()==MUL)
      _pStack.push(L*R);
    else if(_opStack.top()==POW)
      _pStack.push(pow(L,R));
    else return false;
    _opStack.pop();
    return true;
  }
  //data
  std::stack<TOKEN> _opStack;
  std::stack<POLY> _pStack;
};
