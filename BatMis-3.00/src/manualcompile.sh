make clean
make
cd BatMis-3.00/src
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT adler32.o -MD -MP -MF .deps/adler32.Tpo -c -o adler32.o adler32.c -lm
mv -f .deps/adler32.Tpo .deps/adler32.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT compress.o -MD -MP -MF .deps/compress.Tpo -c -o compress.o compress.c -lm
mv -f .deps/compress.Tpo .deps/compress.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT crc32.o -MD -MP -MF .deps/crc32.Tpo -c -o crc32.o crc32.c -lm
mv -f .deps/crc32.Tpo .deps/crc32.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT gzio.o -MD -MP -MF .deps/gzio.Tpo -c -o gzio.o gzio.c -lm
mv -f .deps/gzio.Tpo .deps/gzio.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT uncompr.o -MD -MP -MF .deps/uncompr.Tpo -c -o uncompr.o uncompr.c -lm
mv -f .deps/uncompr.Tpo .deps/uncompr.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT deflate.o -MD -MP -MF .deps/deflate.Tpo -c -o deflate.o deflate.c -lm
mv -f .deps/deflate.Tpo .deps/deflate.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT trees.o -MD -MP -MF .deps/trees.Tpo -c -o trees.o trees.c -lm
mv -f .deps/trees.Tpo .deps/trees.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT zutil.o -MD -MP -MF .deps/zutil.Tpo -c -o zutil.o zutil.c -lm
mv -f .deps/zutil.Tpo .deps/zutil.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT inflate.o -MD -MP -MF .deps/inflate.Tpo -c -o inflate.o inflate.c -lm
mv -f .deps/inflate.Tpo .deps/inflate.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT infback.o -MD -MP -MF .deps/infback.Tpo -c -o infback.o infback.c -lm
mv -f .deps/infback.Tpo .deps/infback.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT inftrees.o -MD -MP -MF .deps/inftrees.Tpo -c -o inftrees.o inftrees.c -lm
mv -f .deps/inftrees.Tpo .deps/inftrees.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT inffast.o -MD -MP -MF .deps/inffast.Tpo -c -o inffast.o inffast.c -lm
mv -f .deps/inffast.Tpo .deps/inffast.Po
rm -f libz.a
ar cru libz.a adler32.o compress.o crc32.o gzio.o uncompr.o deflate.o trees.o zutil.o inflate.o infback.o inftrees.o inffast.o 
ranlib libz.a
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT bwtformatdb.o -MD -MP -MF .deps/bwtformatdb.Tpo -c -o bwtformatdb.o bwtformatdb.c -lm
mv -f .deps/bwtformatdb.Tpo .deps/bwtformatdb.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT BWT.o -MD -MP -MF .deps/BWT.Tpo -c -o BWT.o BWT.c -lm
mv -f .deps/BWT.Tpo .deps/BWT.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT BWTConstruct.o -MD -MP -MF .deps/BWTConstruct.Tpo -c -o BWTConstruct.o BWTConstruct.c -lm
mv -f .deps/BWTConstruct.Tpo .deps/BWTConstruct.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT MiscUtilities.o -MD -MP -MF .deps/MiscUtilities.Tpo -c -o MiscUtilities.o MiscUtilities.c -lm
mv -f .deps/MiscUtilities.Tpo .deps/MiscUtilities.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT MemManager.o -MD -MP -MF .deps/MemManager.Tpo -c -o MemManager.o MemManager.c -lm
mv -f .deps/MemManager.Tpo .deps/MemManager.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT TextConverter.o -MD -MP -MF .deps/TextConverter.Tpo -c -o TextConverter.o TextConverter.c -lm
mv -f .deps/TextConverter.Tpo .deps/TextConverter.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT r250.o -MD -MP -MF .deps/r250.Tpo -c -o r250.o r250.c -lm
mv -f .deps/r250.Tpo .deps/r250.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT QSufSort.o -MD -MP -MF .deps/QSufSort.Tpo -c -o QSufSort.o QSufSort.c -lm
mv -f .deps/QSufSort.Tpo .deps/QSufSort.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT iniparser.o -MD -MP -MF .deps/iniparser.Tpo -c -o iniparser.o iniparser.c -lm
mv -f .deps/iniparser.Tpo .deps/iniparser.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT inistrlib.o -MD -MP -MF .deps/inistrlib.Tpo -c -o inistrlib.o inistrlib.c -lm
mv -f .deps/inistrlib.Tpo .deps/inistrlib.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT dictionary.o -MD -MP -MF .deps/dictionary.Tpo -c -o dictionary.o dictionary.c -lm
mv -f .deps/dictionary.Tpo .deps/dictionary.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT DNACount.o -MD -MP -MF .deps/DNACount.Tpo -c -o DNACount.o DNACount.c -lm
mv -f .deps/DNACount.Tpo .deps/DNACount.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT Timing.o -MD -MP -MF .deps/Timing.Tpo -c -o Timing.o Timing.c -lm
mv -f .deps/Timing.Tpo .deps/Timing.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT Socket.o -MD -MP -MF .deps/Socket.Tpo -c -o Socket.o Socket.c -lm
mv -f .deps/Socket.Tpo .deps/Socket.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT HSP.o -MD -MP -MF .deps/HSP.Tpo -c -o HSP.o HSP.c -lm
mv -f .deps/HSP.Tpo .deps/HSP.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT HSPstatistic.o -MD -MP -MF .deps/HSPstatistic.Tpo -c -o HSPstatistic.o HSPstatistic.c -lm
mv -f .deps/HSPstatistic.Tpo .deps/HSPstatistic.Po
gcc -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2 -MT karlin.o -MD -MP -MF .deps/karlin.Tpo -c -o karlin.o karlin.c -lm
mv -f .deps/karlin.Tpo .deps/karlin.Po
gcc -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2  -g -O2   -o bwtformatdb bwtformatdb.o BWT.o BWTConstruct.o MiscUtilities.o MemManager.o TextConverter.o r250.o QSufSort.o iniparser.o inistrlib.o dictionary.o DNACount.o Timing.o Socket.o HSP.o HSPstatistic.o karlin.o   -lm
g++ -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2   -g -O2 -MT reverse.o -MD -MP -MF .deps/reverse.Tpo -c -o reverse.o reverse.cpp -lm
mv -f .deps/reverse.Tpo .deps/reverse.Po
g++ -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2   -g -O2   -o reverse reverse.o   -lm
g++ -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2   -g -O2 -MT batman.o -MD -MP -MF .deps/batman.Tpo -c -o batman.o batman.cpp -lm
mv -f .deps/batman.Tpo .deps/batman.Po
g++ -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2   -g -O2   -o batman batman.o BWT.o BWTConstruct.o MiscUtilities.o MemManager.o TextConverter.o r250.o QSufSort.o iniparser.o inistrlib.o dictionary.o DNACount.o Timing.o Socket.o HSP.o HSPstatistic.o karlin.o libz.a  -lm
g++ -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2   -g -O2 -MT batdecode.o -MD -MP -MF .deps/batdecode.Tpo -c -o batdecode.o batdecode.cpp -lm
mv -f .deps/batdecode.Tpo .deps/batdecode.Po
g++ -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2   -g -O2   -o batdecode batdecode.o BWT.o BWTConstruct.o MiscUtilities.o MemManager.o TextConverter.o r250.o QSufSort.o iniparser.o inistrlib.o dictionary.o DNACount.o Timing.o Socket.o HSP.o HSPstatistic.o karlin.o libz.a  -lm
g++ -DHAVE_CONFIG_H -I. -I..    -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2   -g -O2 -MT filter.o -MD -MP -MF .deps/filter.Tpo -c -o filter.o filter.cpp -lm
mv -f .deps/filter.Tpo .deps/filter.Po
g++ -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2   -g -O2   -o filter filter.o   -lm
make copy
