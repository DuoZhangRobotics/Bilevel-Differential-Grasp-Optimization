############################################################################
## So@GUI@/data/bin/Makefile.am
## $Id$

BUILT_SOURCES = iv2h.sh

EXTRA_DIST = \
	data2iv.sh \
	iv2h.sh.in \
	material-node.sh

all-local: iv2h.sh

DISTCLEANFILES = $(BUILT_SOURCES)

############################################################################

iv2h.sh: $(top_srcdir)/bin/iv2h.sh.in
	@( cd $(top_builddir); $(MAKE) built-sources )

############################################################################
