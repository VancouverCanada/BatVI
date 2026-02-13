EXTRA_DIST =  batman.ini 
AM_CFLAGS = -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm
AM_CXXFLAGS = -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm 
RANLIB=ranlib
#AM_CPPFLAGS = -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -lm

noinst_LIBRARIES = libz.a
libz_a_SOURCES = adler32.c compress.c crc32.c gzio.c uncompr.c deflate.c trees.c \
       zutil.c inflate.c infback.c inftrees.c inffast.c\
       zlib.h zconf.h crc32.h  deflate.h zutil.h inftrees.h inflate.h inffast.h trees.h inffixed.h

bin_PROGRAMS = bwtformatdb reverse batman batdecode filter
bwtformatdb_SOURCES= bwtformatdb.c BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c\
 BWT.h TypeNLimit.h MemManager.h TextConverter.h HSP.h MiscUtilities.h DNACount.h r250.h HSPstatistic.h\
 BWTConstruct.h QSufSort.h r250.h dictionary.h inistrlib.h iniparser.h Timing.h Socket.h karlin.h
batman_SOURCES= batman.cpp BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c 
batman_LDADD = libz.a
batdecode_SOURCES= batdecode.cpp BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c r250.c QSufSort.c\
 iniparser.c inistrlib.c dictionary.c DNACount.c Timing.c Socket.c HSP.c HSPstatistic.c karlin.c
batdecode_LDADD = libz.a
reverse_SOURCES = reverse.cpp 
filter_SOURCES = filter.cpp 
copy:
	if [ -d "../bin" ]; then echo ../bin exists; else mkdir ../bin; fi
	cp bwtformatdb ../bin
	cp reverse ../bin
	cp batman ../bin
	cp batdecode ../bin
	cp filter ../bin
