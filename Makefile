JVMFLAGS = -cp build
JCFLAGS = -cp src -d build
JC = javac
JVM= java 
CLASSES=$(wildcard *.java)
PKG=mnkgame
MAIN=$(PKG).MNKGame 
SIZE=3

.SUFFIXES: .java .class
random: build
	$(JVM) $(JVMFLAGS) $(MAIN) $(SIZE) $(SIZE) $(SIZE) $(PKG).RandomPlayer

build: $(CLASSES)
	mkdir build
	$(JC) $(JCFLAGS) src/mnkgame/*.java

clean:
	rm -rf build
