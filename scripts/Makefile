../src/sequences/alphabet.hpp:alphabet.py
	${eval TMPFILE := ${shell mktemp /tmp/alphabet.hpp}}
	./alphabet.py > ${TMPFILE}
	clang-format ${TMPFILE} > $@
