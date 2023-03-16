bin_PROGRAMS = kh

dSFMT = ./head/dSFMT-common.h \
				./head/dSFMT-params.h \
				./head/dSFMT-params521.h \
				./head/dSFMT-params1279.h \
				./head/dSFMT-params2203.h \
				./head/dSFMT-params4253.h \
				./head/dSFMT-params11213.h \
				./head/dSFMT-params19937.h \
				./head/dSFMT-params44497.h \
				./head/dSFMT-params86243.h \
				./head/dSFMT-params132049.h \
				./head/dSFMT-params216091.h \
				./head/dSFMT.h \
				./lib/dSFMT.c

libs = 	./lib/spin_manip.c \
				./lib/in_out.c \
				./lib/initialize.c \
				./lib/measures.c \
				./lib/updates.c \
				./lib/utility.c

kh_SOURCES = 	./src/main.c \
							./head/head_&_structures.h \
							$(dSFMT) \
							$(libs)

CC = gcc
AM_CFLAGS = -O3 -msse2 -DDSFMT_MEXP=19937 -DHAVE_SSE2 -fPIE
