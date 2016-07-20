AM_YFLAGS = -d
BUILT_SOURCES = newick_parser.h

noinst_LTLIBRARIES = libnetabc.la
libnetabc_la_SOURCES = util.h util.c stats.c stats.h smc.c smc.h tree.h tree.c treestats.h treestats.c simulate.h simulate.c mmpp.h mmpp.c newick_parser.h newick_parser.y newick_lexer.l $(top_builddir)/igraph/include/igraph.h
libnetabc_la_CFLAGS = $(PTHREAD_CFLAGS) $(GSL_CFLAGS) $(LIBR_CFLAGS) $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include -I$(top_srcdir)/c-cmaes
libnetabc_la_LIBADD = $(PTHREAD_LIBS) $(GSL_LIBS) $(LIBR_LIBS) $(WARN_LDFLAGS) $(top_builddir)/igraph/src/libigraph.la $(top_builddir)/c-cmaes/libcmaes.la -lm -lstdc++ -lgmp -lxml2

bin_PROGRAMS = nettree treekernel netabc treestat pcbr

netabc_SOURCES = netabc.c
netabc_CFLAGS = $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include
netabc_LDADD = $(WARN_LDFLAGS) libnetabc.la

nettree_SOURCES = nettree.c 
nettree_CFLAGS = $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include
nettree_LDADD = $(WARN_LDFLAGS) libnetabc.la

treekernel_sources = treekernel.c
treekernel_CFLAGS = $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include
treekernel_LDADD = $(WARN_LDFLAGS) libnetabc.la

treestat_sources = treestat.c
treestat_CFLAGS = $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include
treestat_LDADD = $(WARN_LDFLAGS) libnetabc.la

pcbr_sources = pcbr.c
pcbr_CFLAGS = $(WARN_CFLAGS) -I$(top_builddir)/igraph/include -I$(top_srcdir)/igraph/include -I$(top_srcdir)/c-cmaes
pcbr_LDADD = $(WARN_LDFLAGS) libnetabc.la
