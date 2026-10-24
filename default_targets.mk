include ${GTTL}/compile_targets.mk

.PHONY:code_check
code_check:code_check_cpp code_check_py

# as code_check.py requires at least one argument, we only call it if there is
# a least one file with the corresponding suffixes in the current directory
.PHONY:code_check_cpp
code_check_cpp:
	@cpp_files=$$(find . -name '*.[hc]pp') && test -z "$${cpp_files}" || ${GTTL}/scripts/code_check.py -wt $${cpp_files}

.PHONY:code_check_py
code_check_py:
	@py_files=$$(find . -name '*.py') && test -z "$${py_files}" || ${GTTL}/scripts/code_check.py -wt $${py_files}

.PHONY:clean
clean:
	@${RM} -r *.[oxd] */*.[oxd] tmp.* TMP* *.fast[aq] *.fast[aq].gz
	@${RM} -r sfx.* *.x.dSYM/ __pycache__ simd32_seqs
# The following are test-files generated in tools/ntcard
	@${RM} eco29.fasta soil_unique.fasta
	@${RM} 70x_161nt_phread64.fastq.gz SRR19536726_1_1000.fastq
# The following is for tools/suffixarrays
	@${GTTL}/scripts/cleanpp.sh
	@${RM} -r sa_induced.x.dSYM

.PHONY:cxx_options
cxx_options:
	@echo $(CXX) $(TIME_OPTION) -DCXX_FLAGS=\"'$(CXXFLAGS) $(CFLAGS) $(CPPFLAGS)'\" -c "<file.cpp>" -o "<file.o>" $(CXXFLAGS) $(CFLAGS) $(CPPFLAGS) -MT "<file.o>" -MMD -MP -MF "<file.d>"
	@echo $(LD) $(LDFLAGS) "<file.o>" -o "<file.x>" $(LDLIBS)
