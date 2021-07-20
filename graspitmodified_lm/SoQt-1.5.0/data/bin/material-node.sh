#!/bin/sh
############################################################################
# material-node.sh
#
# This script takes the environment variables ambient{0,1,2},
# diffuse{0,1,2}, specular{0,1,2}, emissive{0,1,2}, shininess0, and
# transparency0, and outputs an Inventor file containing a Material
# node with the given values.  Unspecified values default to 0.
#
# TODO
# - [nit] get final trailing 0 removed from floats
#
# Authors:
#   Lars J. Aas <larsa@sim.no>
#

if test $# -ne 0; then
  echo "Usage: $programname"
  exit 1
fi

# set up default values
: ${ambient0=0} ${ambient1=0} ${ambient2=0}
: ${diffuse0=0} ${diffuse1=0} ${diffuse2=0}
: ${specular0=0} ${specular1=0} ${specular2=0}
: ${emissive0=0} ${emissive1=0} ${emissive2=0}
: ${shininess0=0}
: ${transparency0=0}

cat $0 | sed -e '1,/^DATA$/ d' \
             -e "s/\<ambient0\>/${ambient0}/" \
             -e "s/\<ambient1\>/${ambient1}/" \
             -e "s/\<ambient2\>/${ambient2}/" \
             -e "s/\<diffuse0\>/${diffuse0}/" \
             -e "s/\<diffuse1\>/${diffuse1}/" \
             -e "s/\<diffuse2\>/${diffuse2}/" \
             -e "s/\<specular0\>/${specular0}/" \
             -e "s/\<specular1\>/${specular1}/" \
             -e "s/\<specular2\>/${specular2}/" \
             -e "s/\<emissive0\>/${emissive0}/" \
             -e "s/\<emissive1\>/${emissive1}/" \
             -e "s/\<emissive2\>/${emissive2}/" \
             -e "s/\<shininess0\>/${shininess0}/" \
             -e "s/\<transparency0\>/${transparency0}/" |
         sed -e 's/\([^0-9]\)\./\10./g' -e 's/00*\>/0/g'

exit 0

DATA
#Inventor V1.0 ascii

Material {
	ambientColor	ambient0 ambient1 ambient2
	diffuseColor	diffuse0 diffuse1 diffuse2
	specularColor	specular0 specular1 specular2
	emissiveColor	emissive0 emissive1 emissive2
	shininess	shininess0
	transparency	transparency0
}
