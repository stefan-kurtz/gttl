TOOL_DIRS=$(sort $(dir $(wildcard tools/*/)))

all:
	@make -C testsuite test
	$(foreach tool,$(TOOL_DIRS), @(make -C $(tool) test))

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
