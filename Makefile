all:
	make -C testsuite test

debug:
	make -C testsuite debug=yes test

.PHONY:clean	
clean:
	make -C testsuite clean

.PHONY:tags
tags:
	./scripts/ctags.sh

.PHONY:code_check
code_check:
	code_check.py -wt `find . -name '*.[ch]pp'`
