gVMFLAGS = -cp build
JCFLAGS = -cp src -d build
JC = javac
JVM= java 
JVMFLAGS = -cp build
CLASSES=$(wildcard *.java)

PKG=mnkgame
MAIN=$(PKG).MNKGame 
PLAYERS=$(PKG).players
SIZE:=3

.SUFFIXES: .java .class
.PHONY: build

albert: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).AlbertEinstein

random: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).QuasiRandomPlayer

build: $(CLASSES)
	mkdir -p build
	$(JC) $(JCFLAGS) src/mnkgame/*.java src/mnkgame/players/*.java

clean:
	rm -rf build
