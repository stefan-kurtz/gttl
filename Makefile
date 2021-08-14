.PHONY:code_check
code_check:
	code_check.py -wt `find . -name '*.[ch]pp'`
