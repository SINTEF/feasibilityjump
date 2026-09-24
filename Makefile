CXXFLAGS ?= -O3 -Wall
override CXXFLAGS += -std=c++11
CPPFLAGS += -I"$(XPRESSDIR)/include"
LDFLAGS  += -L"$(XPRESSDIR)/lib"
LDLIBS   += -lxprs -lpthread

all: xpress_fj

xpress_fj: xpress_fj.cc feasibilityjump.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)
