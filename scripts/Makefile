TARGET_DIR=../src/sequences
.PHONY:../src/sequences/alphabet.hpp
${TARGET_DIR}/alphabet.hpp:
	${eval TMPFILE := ${shell mktemp /tmp/alphabet.hpp.XXXXXX}}
	./alphabet.py > ${TMPFILE}
	cd ${TARGET_DIR} && clang-format --assume-filename=alphabet.cpp ${TMPFILE} > alphabet.hpp
	rm -f ${TMPFILE}
