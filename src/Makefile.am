AM_YFLAGS = -d
BUILT_SOURCES = newick_parser.h

lib_LTLIBRARIES = libnetabc.la
libnetabc_la_SOURCES = util.h util.c stats.c stats.h smc.c smc.h tree.h tree.c treestats.h treestats.c simulate.h simulate.c mmpp.h mmpp.c newick_parser.h newick_parser.y newick_lexer.l $(top_builddir)/igraph/include/igraph.h
libnetabc_la_CFLAGS = $(PTHREAD_CFLAGS) $(GSL_CFLAGS) $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include -I$(top_srcdir)/c-cmaes
libnetabc_la_LIBADD = $(WARN_LDFLAGS) $(top_builddir)/igraph/src/libigraph.la $(top_builddir)/c-cmaes/libcmaes.la $(PTHREAD_LIBS) $(GSL_LIBS) -lm -lstdc++ -lgmp -lxml2

include_HEADERS = netabc.h
pkginclude_HEADERS = smc.h util.h stats.h tree.h treestats.h simulate.h mmpp.h newick_parser.h

bin_PROGRAMS = nettree treekernel netabc treestat pcbr

netabc_SOURCES = netabc.c
netabc_CFLAGS = $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include $(GSL_CFLAGS)
netabc_LDADD = $(WARN_LDFLAGS) libnetabc.la $(GSL_LIBS)

nettree_SOURCES = nettree.c 
nettree_CFLAGS = $(GSL_CFLAGS) $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include
nettree_LDADD = libnetabc.la $(GSL_LIBS) $(WARN_LDFLAGS)

treekernel_sources = treekernel.c
treekernel_CFLAGS = $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include
treekernel_LDADD = $(WARN_LDFLAGS) libnetabc.la

treestat_sources = treestat.c
treestat_CFLAGS = $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include $(GSL_CFLAGS)
treestat_LDADD = $(WARN_LDFLAGS) libnetabc.la $(GSL_LIBS)

pcbr_sources = pcbr.c
pcbr_CFLAGS = $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include -I$(top_srcdir)/c-cmaes
pcbr_LDADD = $(WARN_LDFLAGS) libnetabc.la
