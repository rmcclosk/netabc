AM_YFLAGS = -d
BUILT_SOURCES = newick_parser.h

noinst_LTLIBRARIES = libnetabc.la
libnetabc_la_SOURCES = util.h util.c stats.c stats.h smc.c smc.h tree.h tree.c treestats.h treestats.c simulate.h simulate.c newick_parser.h newick_parser.y newick_lexer.l
libnetabc_la_CFLAGS = $(PTHREAD_CFLAGS) $(GSL_CFLAGS) $(LIBR_CFLAGS) $(WARN_CFLAGS)
libnetabc_la_LIBADD = $(PTHREAD_LIBS) $(GSL_LIBS) $(LIBR_LIBS) $(WARN_LDFLAGS)

bin_PROGRAMS = nettree treekernel netabc treestat

netabc_SOURCES = netabc.c
netabc_CFLAGS = $(WARN_CFLAGS)
netabc_LDADD = $(WARN_LDFLAGS) libnetabc.la

nettree_SOURCES = nettree.c 
nettree_CFLAGS = $(WARN_CFLAGS) 
nettree_LDADD = $(WARN_LDFLAGS) libnetabc.la

treekernel_sources = treekernel.c
treekernel_CFLAGS = $(WARN_CFLAGS) 
treekernel_LDADD = $(WARN_LDFLAGS) libnetabc.la

treestat_sources = treestat.c
treestat_CFLAGS = $(WARN_CFLAGS) 
treestat_LDADD = $(WARN_LDFLAGS) libnetabc.la
