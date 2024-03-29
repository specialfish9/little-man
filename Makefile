JVMFLAGS = -cp build
JCFLAGS = -cp src -d build
JC = javac
JVM= java 
JVMFLAGS = -cp build
CLASSES=$(wildcard *.java)

PKG=mnkgame
MAIN=$(PKG).MNKGame 
PLAYERS=$(PKG).players
M:=4
N:=4
K:=3
REPS:=5
TIME:=10
BEST:=LittleMan
OLD:=LittleBoy

.SUFFIXES: .java .class
.PHONY: build

best: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).$(BEST)

human: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K)

vs: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).$(BEST) $(PLAYERS).$(OLD)

vs1: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(M) $(N) $(K) $(PLAYERS).$(OLD) $(PLAYERS).$(BEST)

test: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(M) $(N) $(K) $(PLAYERS).$(BEST) $(PLAYERS).$(OLD) -v -r $(REPS)

ttest: build
	$(JVM) $(JVMFLAGS) mnkgame.MNKPlayerTester $(M) $(N) $(K) $(PLAYERS).$(BEST) $(PLAYERS).$(BEST) -v -r $(REPS)

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
