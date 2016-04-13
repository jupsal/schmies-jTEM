SRCDIRS=../src

BINDIR=../classes

JUNIT=$(wildcard junit*.jar)

LIBDIR=../lib

JAVACOPTS=-target 1.6 -source 1.6

SOURCEFILES=$(shell find $(SRCDIRS) -name '*.java')

DEPNAMES=$(shell cat ../dependencies.txt 2> /dev/null | grep -v '^\#' )
DEPS=$(patsubst %,$(LIBDIR)/%.jar, $(DEPNAMES))

CLASSES = \
		  justASchottkyGroup.java \
		  SchottkyGroupExample.java \
		  GetSchottkyData.java \
		  EmployeeTest.java

.PHONY: jeremy
jeremy: $(CLASSES:.java=.class)
$(CLASSES:.java=.class): $(CLASSES:.java=.class)
	javac $(JAVACOPTS) \
		-classpath ../lib/mfc.jar:../lib/blas.jar:../lib/numericalMethods.jar:../classes/ \
		-sourcepath ../src/ \
		-d . \
	   	$(CLASSES)
	
	
clean: ; rm -rf *.class
