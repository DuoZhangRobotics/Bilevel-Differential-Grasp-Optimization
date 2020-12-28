#ifndef LEVI_CIVITA_H
#define LEVI_CIVITA_H

#include <CommonFile/Config.h>

PRJ_BEGIN

template <sizeType i,sizeType ROrC>struct LeviCivitaPositive;
template <sizeType i,sizeType ROrC>struct LeviCivitaNegative;
//i=0
template <>struct LeviCivitaPositive<0,0>{static const sizeType value=1;};
template <>struct LeviCivitaPositive<0,1>{static const sizeType value=2;};
template <>struct LeviCivitaNegative<0,0>{static const sizeType value=2;};
template <>struct LeviCivitaNegative<0,1>{static const sizeType value=1;};
//i=1
template <>struct LeviCivitaPositive<1,0>{static const sizeType value=2;};
template <>struct LeviCivitaPositive<1,1>{static const sizeType value=0;};
template <>struct LeviCivitaNegative<1,0>{static const sizeType value=0;};
template <>struct LeviCivitaNegative<1,1>{static const sizeType value=2;};
//i=2
template <>struct LeviCivitaPositive<2,0>{static const sizeType value=0;};
template <>struct LeviCivitaPositive<2,1>{static const sizeType value=1;};
template <>struct LeviCivitaNegative<2,0>{static const sizeType value=1;};
template <>struct LeviCivitaNegative<2,1>{static const sizeType value=0;};

PRJ_END

#endif
