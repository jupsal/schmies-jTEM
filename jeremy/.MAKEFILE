SRCDIRS=../src

BINDIR=classes

SOURCEFILES=$(shell find $(SRCDIRS) -name '*.java')

DEPNAMES=$(shell cat ../dependencies.txt 2> /dev/null | grep -v '^\#' )
DEPS=$(patsubst %,$(LIBDIR)/%.jar, $(DEPNAMES))


.PHONY: jeremy
jeremy: justASchottkyGroup.class
justASchottkyGroup.class: justASchottkyGroup.java ; javac -d . -classpath . justASchottkyGroup.java

clean: ; rm -f *.class
