#!/bin/sh
############################################################################
# data2iv.sh
#
# This script generates a set of material node Inventor files.
# The input file is named <class>.dat and contains a set of lines following
# the following format:
#
# ambient:diffuse:emissive:specular:shininess:transparency
#
# The four first fields are color values, and must be specified either on
# the "#RRGGBB" format or the more precise "#RRRRGGGGBBBB" format.
# Hexadecimal of course.
#
# The last two fields are single value fields, either specified through hex
# "#XX" (scaled by 0xff) or "#XXXX" (scaled by 0xffff), or as a decimal
# number (also scaled by 65535) or as an unscaled floating point number.
#
# Fields can be omitted (they will then default to contain zero-values), as
# long as the field-separating ":"s are still in place.
#
# For something cool, check out how hex-values are translated to decimal
# numbers ;)
#
# TODO:
# - use something other than dc to calculate float values?
#
# Authors:
#   Lars J. Aas <larsa@sim.no>
#

if test $# -ne 1; then
  echo "Usage: $programname <materials-file.dat>"
  exit 1
fi

############################################################################

dehexXX='sed -e "s/\\(.\\)\\(.\\)/\\1 * 16 + \\2/" -e "s/\\([abcdefABCDEF]\\)/1\\1/g" -e "y/abcdefABCDEF/012345012345/"'
dehexXXXX='sed -e "s/\\(.\\)\\(.\\)\\(.\\)\\(.\\)/\\1 * 4096 + \\2 * 256 + \\3 * 16 + \\4/" -e "s/\\([abcdefABCDEF]\\)/1\\1/g" -e "y/abcdefABCDEF/012345012345/"'

precision=8
file=$1
style=`basename $1 .dat`
dir="`echo $0 | sed -e 's/[^\/]*$//g'`"

test ! -e $style && mkdir $style
test ! -e $dir/material-node.sh && exit 1

line=1
lines=`cat $file | wc -l | xargs expr 0 +`
while test $line -le $lines; do
  data=`cat $file | sed -n -e "$line p"`
  save_IFS=$IFS
  IFS=:
  set $data
  IFS=$save_IFS

  ambient=$1
  diffuse=$2
  emissive=$3
  specular=$4
  shininess0=$5
  transparency0=$6

# ambient color
  if test x"$ambient" = "x"; then
    ambient0=0
    ambient1=0
    ambient2=0
  elif test x"`echo $ambient | egrep "^#"`" != x; then
    if test "`echo $ambient | wc -c`" -eq 14; then
      # assuming "#RRRRGGGGBBBB" format
      ambient0=`echo $ambient | sed -e 's/^#\(....\).*$/\1/'`
      ambient0=`echo $ambient0 | eval "$dehexXXXX" | xargs expr`
      ambient0=`echo "$precision k $ambient0 65535 / n" | dc`
      ambient1=`echo $ambient | sed -e 's/^#....\(....\).*$/\1/'`
      ambient1=`echo $ambient1 | eval "$dehexXXXX" | xargs expr`
      ambient1=`echo "$precision k $ambient1 65535 / n" | dc`
      ambient2=`echo $ambient | sed -e 's/^#........\(....\)$/\1/'`
      ambient2=`echo $ambient2 | eval "$dehexXXXX" | xargs expr`
      ambient2=`echo "$precision k $ambient2 65535 / n" | dc`
    elif test "`echo $ambient | wc -c`" -eq 8; then
      # assuming "#RRGGBB" format
      ambient0=`echo $ambient | sed -e 's/^#\(..\).*$/\1/'`
      ambient0=`echo $ambient0 | eval "$dehexXX" | xargs expr`
      ambient0=`echo "$precision k $ambient0 255 / n" | dc`
      ambient1=`echo $ambient | sed -e 's/^#..\(..\).*$/\1/'`
      ambient1=`echo $ambient1 | eval "$dehexXX" | xargs expr`
      ambient1=`echo "$precision k $ambient1 255 / n" | dc`
      ambient2=`echo $ambient | sed -e 's/^#....\(..\)$/\1/'`
      ambient2=`echo $ambient2 | eval "$dehexXX" | xargs expr`
      ambient2=`echo "$precision k $ambient2 255 / n" | dc`
    else
      echo "Error: line $line: format of ambient color field is not supported"
      exit 1
    fi
  else
    echo "Error: line $line: format of ambient color field is not supported"
    exit 1
  fi
#  echo "ambient: '$ambient0' '$ambient1' '$ambient2'"

  # diffuse color
  if test x"$diffuse" = "x"; then
    diffuse0=0
    diffuse1=0
    diffuse2=0
  elif test x"`echo $diffuse | egrep "^#"`" != x; then
    if test "`echo $diffuse | wc -c`" -eq 14; then
      # assuming "#RRRRGGGGBBBB" format
      diffuse0=`echo $diffuse | sed -e 's/^#\(....\).*$/\1/'`
      diffuse0=`echo $diffuse0 | eval "$dehexXXXX" | xargs expr`
      diffuse0=`echo "$precision k $diffuse0 65535 / n" | dc`
      diffuse1=`echo $diffuse | sed -e 's/^#....\(....\).*$/\1/'`
      diffuse1=`echo $diffuse1 | eval "$dehexXXXX" | xargs expr`
      diffuse1=`echo "$precision k $diffuse1 65535 / n" | dc`
      diffuse2=`echo $diffuse | sed -e 's/^#........\(....\)$/\1/'`
      diffuse2=`echo $diffuse2 | eval "$dehexXXXX" | xargs expr`
      diffuse2=`echo "$precision k $diffuse2 65535 / n" | dc`
    elif test "`echo $diffuse | wc -c`" -eq 8; then
      # assuming "#RRGGBB" format
      diffuse0=`echo $diffuse | sed -e 's/^#\(..\).*$/\1/'`
      diffuse0=`echo $diffuse0 | eval "$dehexXX" | xargs expr`
      diffuse0=`echo "$precision k $diffuse0 255 / n" | dc`
      diffuse1=`echo $diffuse | sed -e 's/^#..\(..\).*$/\1/'`
      diffuse1=`echo $diffuse1 | eval "$dehexXX" | xargs expr`
      diffuse1=`echo "$precision k $diffuse1 255 / n" | dc`
      diffuse2=`echo $diffuse | sed -e 's/^#....\(..\)$/\1/'`
      diffuse2=`echo $diffuse2 | eval "$dehexXX" | xargs expr`
      diffuse2=`echo "$precision k $diffuse2 255 / n" | dc`
    else
      echo "Error: line $line: format of diffuse color field is not supported"
      exit 1
    fi
  else
    echo "Error: line $line: format of diffuse color field is not supported"
    exit 1
  fi
#  echo "diffuse: '$diffuse0' '$diffuse1' '$diffuse2'"
  
  if test x"$specular" = "x"; then
    specular0=0
    specular1=0
    specular2=0
  elif test x"`echo $specular | egrep "^#"`" != x; then
    if test "`echo $specular | wc -c`" -eq 14; then
      # assuming "#RRRRGGGGBBBB" format
      specular0=`echo $specular | sed -e 's/^#\(....\).*$/\1/'`
      specular0=`echo $specular0 | eval "$dehexXXXX" | xargs expr`
      specular0=`echo "$precision k $specular0 65535 / n" | dc`
      specular1=`echo $specular | sed -e 's/^#....\(....\).*$/\1/'`
      specular1=`echo $specular1 | eval "$dehexXXXX" | xargs expr`
      specular1=`echo "$precision k $specular1 65535 / n" | dc`
      specular2=`echo $specular | sed -e 's/^#........\(....\)$/\1/'`
      specular2=`echo $specular2 | eval "$dehexXXXX" | xargs expr`
      specular2=`echo "$precision k $specular2 65535 / n" | dc`
    elif test "`echo $specular | wc -c`" -eq 8; then
      # assuming "#RRGGBB" format
      specular0=`echo $specular | sed -e 's/^#\(..\).*$/\1/'`
      specular0=`echo $specular0 | eval "$dehexXX" | xargs expr`
      specular0=`echo "$precision k $specular0 255 / n" | dc`
      specular1=`echo $specular | sed -e 's/^#..\(..\).*$/\1/'`
      specular1=`echo $specular1 | eval "$dehexXX" | xargs expr`
      specular1=`echo "$precision k $specular1 255 / n" | dc`
      specular2=`echo $specular | sed -e 's/^#....\(..\)$/\1/'`
      specular2=`echo $specular2 | eval "$dehexXX" | xargs expr`
      specular2=`echo "$precision k $specular2 255 / n" | dc`
    else
      echo "Error: line $line: format of specular color field is not supported"
      exit 1
    fi
  else
    echo "Error: line $line: format of specular color field is not supported"
    exit 1
  fi
#  echo "specular: '$specular0' '$specular1' '$specular2'"

  # emissive color
  if test x"$emissive" = "x"; then
    emissive0=0
    emissive1=0
    emissive2=0
  elif test x"`echo $emissive | egrep "^#"`" != x; then
    if test "`echo $emissive | wc -c`" -eq 14; then
      # assuming "#RRRRGGGGBBBB" format
      emissive0=`echo $emissive | sed -e 's/^#\(....\).*$/\1/'`
      emissive0=`echo $emissive0 | eval "$dehexXXXX" | xargs expr`
      emissive0=`echo "$precision k $emissive0 65535 / n" | dc`
      emissive1=`echo $emissive | sed -e 's/^#....\(....\).*$/\1/'`
      emissive1=`echo $emissive1 | eval "$dehexXXXX" | xargs expr`
      emissive1=`echo "$precision k $emissive1 65535 / n" | dc`
      emissive2=`echo $emissive | sed -e 's/^#........\(....\)$/\1/'`
      emissive2=`echo $emissive2 | eval "$dehexXXXX" | xargs expr`
      emissive2=`echo "$precision k $emissive2 65535 / n" | dc`
    elif test "`echo $emissive | wc -c`" -eq 8; then
      # assuming "#RRGGBB" format
      emissive0=`echo $emissive | sed -e 's/^#\(..\).*$/\1/'`
      emissive0=`echo $emissive0 | eval "$dehexXX" | xargs expr`
      emissive0=`echo "$precision k $emissive0 255 / n" | dc`
      emissive1=`echo $emissive | sed -e 's/^#..\(..\).*$/\1/'`
      emissive1=`echo $emissive1 | eval "$dehexXX" | xargs expr`
      emissive1=`echo "$precision k $emissive1 255 / n" | dc`
      emissive2=`echo $emissive | sed -e 's/^#....\(..\)$/\1/'`
      emissive2=`echo $emissive2 | eval "$dehexXX" | xargs expr`
      emissive2=`echo "$precision k $emissive2 255 / n" | dc`
    else
      echo "Error: line $line: format of emissive color field is not supported"
      exit 1
    fi
  else
    echo "Error: line $line: format of emissive color field is not supported"
    exit 1
  fi
#  echo "emissive: '$emissive0' '$emissive1' '$emissive2'"

  if test x"$shininess0" = "x"; then
    shininess0=0
  elif test x"`echo $shininess0 | egrep "^#"`" != x; then
    if test `echo "$shininess0" | wc -c` -eq 6; then
      # assuming "#XXXX" format
      shininess0=`echo $shininess0 | sed -e 's/^#\(....\)$/\1/'`
      shininess0=`echo $shininess0 | eval "$dehexXXXX" | xargs expr`
      shininess0=`echo "$precision k $shininess0 65535 / n" | dc`
    elif test `echo "$shininess0" | wc -c` -eq 4; then
      # assuming "#XX" format
      shininess0=`echo $shininess0 | sed -e 's/^#\(..\)$/\1/'`
      shininess0=`echo $shininess0 | eval "$dehexXX" | xargs expr`
      shininess0=`echo "$precision k $shininess0 255 / n" | dc`
    else
      echo "Error: line $line: format of shininess field is not supported"
      exit 1
    fi
  elif test `echo "$shininess0" | egrep -c '^[01]?(\.[0-9]+)?$'` -eq 1; then
    # assuming float
    :
  elif test `echo "$shininess0" | egrep -c '^[0-9]+?$'` -eq 1; then
    # assuming 65536-based number
    shininess0=`echo "$precision k $shininess0 65535 / n" | dc`
  else
    echo "Error: line $line: format of shininess field is not supported"
    exit 1
  fi
#  echo "shininess: '$shininess0'"

  if test x"$transparency0" = "x"; then
    transparency0=0
  elif test x"`echo $transparency0 | egrep "^#"`" != x; then
    if test `echo "$transparency0" | wc -c` -eq 6; then
      # assuming "#XXXX" format
      transparency0=`echo $transparency0 | sed -e 's/^#\(....\)$/\1/'`
      transparency0=`echo $transparency0 | eval "$dehexXXXX" | xargs expr`
      transparency0=`echo "$precision k $transparency0 65535 / n" | dc`
    elif test `echo "$transparency0" | wc -c` -eq 4; then
      # assuming "#XX" format
      transparency0=`echo $transparency0 | sed -e 's/^#\(..\)$/\1/'`
      transparency0=`echo $transparency0 | eval "$dehexXX" | xargs expr`
      transparency0=`echo "$precision k $transparency0 255 / n" | dc`
    else
      echo "Error: line $line: format of transparency field is not supported"
      exit 1
    fi
  elif test `echo "$transparency0" | egrep -c '^[01]?(\.[0-9]+)?$'` -eq 1; then
    # assuming float format
    :
  elif test `echo "$transparency0" | egrep -c '^[1-9][0-9]*$'` -eq 1; then
    # assuming 65535-based decimal
    transparency0=`echo "$precision k $transparency0 65535 / n" | dc`
    :
  else
    echo "Error: line $line: format of transparency field is not supported"
    exit 1
  fi
#  echo "transparency: '$transparency0'"

  export ambient0 ambient1 ambient2
  export diffuse0 diffuse1 diffuse2
  export specular0 specular1 specular2
  export emissive0 emissive1 emissive2
  export shininess0 transparency0

  num=`expr $line - 1`
  echo "generating $style/$style.$num"
  $dir/material-node.sh > $style/$style.$num
  line=`expr $line + 1`
done

