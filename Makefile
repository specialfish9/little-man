JVMFLAGS = -cp build
JCFLAGS = -cp src -d build
JC = javac
JVM= java 
JVMFLAGS = -cp build
CLASSES=$(wildcard *.java)

PKG=mnkgame
MAIN=$(PKG).MNKGame 
PLAYERS=$(PKG).players
M:=5
N:=5
K:=5
REPS:=5
TIME:=2
BEST:=GalileoGalilei
OLD:=FrancoRasetti

.SUFFIXES: .java .class
.PHONY: build

best: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).$(BEST)

vs: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).$(BEST) $(PLAYERS).$(OLD)

vs1: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).$(OLD) $(PLAYERS).$(BEST)

test: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(M) $(N) $(K) $(PLAYERS).$(BEST) $(PLAYERS).$(OLD) -v -r $(REPS)

test1: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(M) $(N) $(K) $(PLAYERS).$(OLD) $(PLAYERS).$(BEST) -v -r $(REPS)

constrained: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(M) $(N) $(K) $(PLAYERS).$(OLD) $(PLAYERS).$(BEST) -v -t $(TIME) -r $(REPS)

constrained1: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(M) $(N) $(K) $(PLAYERS).$(BEST) $(PLAYERS).$(OLD) -v -t $(TIME) -r $(REPS)

franco: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).FrancoRasetti

enrico: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).EnricoFermi

david: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).DavidHilbert

charles: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).CharlesDarwin

bertrand: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).BertrandRussell

albert: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).AlbertEinstein

random: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).QuasiRandomPlayer

build: $(CLASSES)
	mkdir -p build
	$(JC) $(JCFLAGS) src/mnkgame/*.java src/mnkgame/players/*.java

clean:
	rm -rf build
