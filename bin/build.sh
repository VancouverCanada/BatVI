gcc -O0 -msse2 -lm -g -c -o ssw.o ssw.c
g++ -O0 -msse2 -lm -g -o cluster_bp ssw.o cluster_bp.cpp
