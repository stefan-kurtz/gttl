all:
	@make -C testsuite test

debug:
	@make -C testsuite debug=yes test

.PHONY:clean	
clean:
	@make -C testsuite clean
	@make -C tools/ntcard clean

.PHONY:tags
tags:
	@./scripts/ctags.sh

.PHONY:code_check
code_check:
	scripts/code_check.py -wt `find . -name '*.[ch]pp' | grep -v merge_sort.hpp`
	scripts/check_ifndef.py `find . -name '*.hpp'`
	scripts/check_no_c_files.sh
	scripts/check_gitattributes.sh
