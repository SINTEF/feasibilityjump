all: xpress_fj

xpress_fj: xpress_fj.cc feasibilityjump.hh
	g++ -O3 -Wall xpress_fj.cc -o xpress_fj -I${XPRESSDIR}/include -L${XPRESSDIR}/lib/ -lxprs -lpthread
