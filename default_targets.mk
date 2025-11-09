include ${GTTL}/compile_targets.mk

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
