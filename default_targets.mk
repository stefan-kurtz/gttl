%.o: %.cpp
	$(CXX) $(TIME_OPTION) -DCXX_FLAGS=\"'$(CXXFLAGS) $(CFLAGS) $(CPPFLAGS)'\" -c $< -o $@ $(CXXFLAGS) $(CFLAGS) $(CPPFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

# In the simplest case, an executable consists of 1 object.
%.x: %.o
	$(LD) $(LDFLAGS) $< -o $@ $(LDLIBS)

.DEFAULT_GOAL := all

.PHONY:tags
tags:
	@ctags -w --c++-kinds=+p --fields=+iaSKlm --extra=+q `ls *.[hc]pp` \
        `find ${GTTL}/src -name '*.hpp'` > $@

.PHONY:code_check
code_check:
	@${GTTL}/scripts/code_check.py -wt `find . -name '*.[hc]pp'`
	@${GTTL}/scripts/code_check.py -wt `find . -name '*.py'`

.PHONY:clean
clean:
	@${RM} -r *.[oxd] */*.[oxd] tmp.* TMP* *.fast[aq] *.fast[aq].gz
	@${RM} -r sfx.* *.x.dSYM/ __pycache__ simd32_seqs
# The following are testing-files generated in tools/ntcard
	@${RM} eco29.fasta soil_unique.fasta 70x_161nt_phread64.fastq.gz SRR19536726_1_1000.fastq
# The following is for tools/suffixarrays
	@${GTTL}/scripts/cleanpp.sh
	@${RM} -r sa_induced.x.dSYM
