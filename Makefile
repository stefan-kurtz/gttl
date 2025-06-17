include config.mk

.PHONY:all
all:
	@make -C testsuite test
	@make -C tools/chaining test
	@make -C tools/ntcard test
	@make -C tools/suffixarrays test
	@make -C tools/swalign test
	@make -C tools/unwords test

.PHONY:debug
debug:
	@make -C testsuite debug=yes test
	@make -C tools/chaining debug=yes test
	@make -C tools/ntcard debug=yes test
	@make -C tools/suffixarrays debug=yes test
	@make -C tools/swalign debug=yes test
	@make -C tools/unwords debug=yes test

.PHONY:clean	
clean:
	@make -C testsuite clean
	@make -C tools/chaining clean
	@make -C tools/ntcard clean
	@make -C tools/suffixarrays clean
	@make -C tools/swalign clean
	@make -C tools/unwords clean

.PHONY:tags
tags:
	@./scripts/ctags.sh

.PHONY:doc
doc:
	@make -C doc

.PHONY:code_check
code_check:
	scripts/code_check.py -wt `find . -name '*.[ch]pp' | grep -v merge_sort.hpp`
	scripts/check_ifndef.py `find . -name '*.hpp'`
	scripts/check_no_c_files.sh
	scripts/check_gitattributes.sh
	scripts/check_portability.sh
	scripts/check_all_sources_tested.py
