JVMFLAGS = -cp build
JCFLAGS = -cp src -d build
JC = javac
JVM= java 
CLASSES=$(wildcard *.java)

PKG=mnkgame
MAIN=$(PKG).MNKGame 
PLAYERS=$(PKG).players
SIZE:=3
PLAYER:=QuasiRandomPlayer

.SUFFIXES: .java .class
.PHONY: build

random: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).$(PLAYER)

build: $(CLASSES)
	mkdir -p build
	$(JC) $(JCFLAGS) src/gametree/*.java src/mnkgame/*.java src/mnkgame/players/*.java

clean:
	rm -rf build
