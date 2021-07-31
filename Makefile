gVMFLAGS = -cp build
JCFLAGS = -cp src -d build
JC = javac
JVM= java 
JVMFLAGS = -cp build
CLASSES=$(wildcard *.java)

PKG=mnkgame
MAIN=$(PKG).MNKGame 
PLAYERS=$(PKG).players
SIZE:=4
K:=3
REPS:=5
BEST:=AlbertEinstein
OLD:=VinDiesel

.SUFFIXES: .java .class
.PHONY: build

best: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(K) $(PLAYERS).$(BEST)

test: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(SIZE) $(SIZE) $(K) $(PLAYERS).$(OLD) $(PLAYERS).$(BEST) -v -r $(REPS)

test1: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(SIZE) $(SIZE) $(K) $(PLAYERS).$(BEST) $(PLAYERS).$(OLD) -v -r $(REPS)

constrained: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(SIZE) $(SIZE) $(K) $(PLAYERS).$(OLD) $(PLAYERS).$(BEST) -v -t 1 -r $(REPS)

constrained1: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(SIZE) $(SIZE) $(K) $(PLAYERS).$(BEST) $(PLAYERS).$(OLD) -v -t 1 -r $(REPS)

albert: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(K) $(PLAYERS).AlbertEinstein

zara: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).Zarathustra

vin: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).VinDiesel

random: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(K) $(PLAYERS).QuasiRandomPlayer

build: $(CLASSES)
	mkdir -p build
	$(JC) $(JCFLAGS) src/mnkgame/*.java src/mnkgame/players/*.java

clean:
	rm -rf build
