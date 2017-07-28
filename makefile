CAPD = /home/marin/Desktop/capd-mp-build/bin/
INCLUDE = `$(CAPD)capd-config --cflags` `$(CAPD)mpcapd-config --cflags` `$(CAPD)capd-gui-config --cflags`
LIBS = `$(CAPD)capd-config --libs` `$(CAPD)mpcapd-config --libs` `$(CAPD)capd-gui-config --libs`

all:
	g++ -O2 -std=c++11 $(INCLUDE) main.cpp $(LIBS) -o main 

