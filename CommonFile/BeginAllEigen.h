//define shorthand names for eigen matrices: templates
#define NAME_EIGEN(type,NAME,size1,size2) \
typedef Eigen::Matrix<type,size1,size2> NAME;

#ifdef NO_CONFLICT
#define NAME_EIGEN_ALLTYPES(NAME,size1,size2) \
NAME_EIGEN(unsigned char,NAME##uc,size1,size2) \
NAME_EIGEN(char,NAME##c,size1,size2) \
NAME_EIGEN(sizeType,NAME##i,size1,size2) \
NAME_EIGEN(scalarF,NAME##f,size1,size2) \
NAME_EIGEN(scalarD,NAME##d,size1,size2)
#else
#define NAME_EIGEN_ALLTYPES(NAME,size1,size2) \
NAME_EIGEN(unsigned char,NAME##uc,size1,size2) \
NAME_EIGEN(char,NAME##c,size1,size2) \
NAME_EIGEN(sizeType,NAME##i,size1,size2) \
NAME_EIGEN(scalarF,NAME##f,size1,size2) \
NAME_EIGEN(scalarD,NAME##d,size1,size2)  \
NAME_EIGEN(scalar,NAME,size1,size2)
#endif

#if defined(FIXED_ONLY)
#define NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,size) \
NAME_EIGEN_ALLTYPES(PREFIX##Row##size,1,size)  \
NAME_EIGEN_ALLTYPES(PREFIX##Vec##size,size,1)
#elif defined(NON_FIXED_ONLY)
#define NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,size) \
NAME_EIGEN_ALLTYPES(PREFIX##Mat##size##X,size,-1)  \
NAME_EIGEN_ALLTYPES(PREFIX##Mat##X##size,-1,size)
#else
#define NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,size) \
NAME_EIGEN_ALLTYPES(PREFIX##Row##size,1,size) \
NAME_EIGEN_ALLTYPES(PREFIX##Vec##size,size,1) \
NAME_EIGEN_ALLTYPES(PREFIX##Mat##size##X,size,-1) \
NAME_EIGEN_ALLTYPES(PREFIX##Mat##X##size,-1,size) \
NAME_EIGEN_ALLTYPES(PREFIX##Mat##size,size,size)
#endif

#if defined(FIXED_ONLY)
#define NAME_EIGEN_ALLTYPES_BASIC(PREFIX)
#else
#define NAME_EIGEN_ALLTYPES_BASIC(PREFIX) \
NAME_EIGEN_ALLTYPES(PREFIX##Mat,-1,-1)  \
NAME_EIGEN_ALLTYPES(PREFIX##Col,-1,1)
#endif

#define NAME_EIGEN_MAT_ALLTYPES(PREFIX,size1,size2)  \
NAME_EIGEN_ALLTYPES(PREFIX##Mat##size1##X##size2,size1,size2)

//define shorthand names for eigen matrices: all types
#define NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE(PREFIX)  \
NAME_EIGEN_ALLTYPES_BASIC(PREFIX) \
NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,18) \
NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,12) \
NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,9)  \
NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,7)  \
NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,6)  \
NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,4)  \
NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,3)  \
NAME_EIGEN_ROWCOL_ALLTYPES(PREFIX,2)  \

#define NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE(PREFIX) \
NAME_EIGEN_MAT_ALLTYPES(PREFIX,9,18)  \
NAME_EIGEN_MAT_ALLTYPES(PREFIX,9,12)  \
NAME_EIGEN_MAT_ALLTYPES(PREFIX,12,12) \
NAME_EIGEN_MAT_ALLTYPES(PREFIX,9,9) \
NAME_EIGEN_MAT_ALLTYPES(PREFIX,6,6) \
NAME_EIGEN_MAT_ALLTYPES(PREFIX,4,4) \
NAME_EIGEN_MAT_ALLTYPES(PREFIX,3,4) \
NAME_EIGEN_MAT_ALLTYPES(PREFIX,3,3) \
NAME_EIGEN_MAT_ALLTYPES(PREFIX,2,2)

#define NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(PREFIX,POSTFIX,type)  \
typedef Eigen::Quaternion<type> PREFIX##Quat##POSTFIX; \
typedef Eigen::Translation<type,3> PREFIX##Trans##POSTFIX; \
typedef Eigen::Transform<type,3,Eigen::Affine> PREFIX##Affine##POSTFIX;
