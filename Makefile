.PHONY:all
all:
	@make -C testsuite test
	@make -C tools/chaining test
	@make -C tools/ntcard test
	@make -C tools/unwords test

.PHONY:build
build:
	@make -C testsuite
	@make -C tools/chaining
	@make -C tools/ntcard
	@make -C tools/unwords

.PHONY:bear
bear:
	@$(MAKE) clean
	@bear -- $(MAKE) $(MAKEGOALS) build
	@sed -i -e '/"-Werror",/d' compile_commands.json

.PHONY:debug
debug:
	@make -C testsuite debug=yes test
	@make -C tools/chaining debug=yes test
	@make -C tools/ntcard debug=yes test
	@make -C tools/unwords debug=yes test

.PHONY:clean	
clean:
	@make -C testsuite clean
	@make -C tools/chaining clean
	@make -C tools/ntcard clean
	@make -C tools/unwords clean

.PHONY:tidy
tidy:
	@find . -name '*.hpp' -o -name '*.cpp' | \
	 xargs -P$(shell nproc) -n1 \
	  clang-tidy --quiet --config-file=.clang-tidy -p . \
	  2>/dev/null

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
