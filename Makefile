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
REPS:=5
BEST:=BertrandRussell
OLD:=AlbertEinstein

.SUFFIXES: .java .class
.PHONY: build

best: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).$(BEST)

test: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).$(OLD) $(PLAYERS).$(BEST) -v -r $(REPS)

test1: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).$(BEST) $(PLAYERS).$(OLD) -v -r $(REPS)

constrained: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).$(OLD) $(PLAYERS).$(BEST) -v -t 1 -r $(REPS)

constrained1: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).$(BEST) $(PLAYERS).$(OLD) -v -t 1 -r $(REPS)

albert: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).AlbertEinstein

random: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(SIZE) $(PLAYERS).QuasiRandomPlayer

build: $(CLASSES)
	mkdir -p build
	$(JC) $(JCFLAGS) src/mnkgame/*.java src/mnkgame/players/*.java

clean:
	rm -rf build
