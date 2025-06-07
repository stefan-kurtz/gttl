all:smrt100-2000.v2.tsv \
    smrt100-2000.v2.best2.tsv \
    smrt100-2000-first.fna \
    smrt100-2000-last.fna \
    query-protein_f5.fsa \
    Q6GZS4.fsa \
    ref-query-protein-v1.tsv

smrt100-2000.v2.tsv:smrt100-2000.a2.tsv
	@cat $< | cut -d '	' -f 1-8 | sed -e 's/, cigar//'  > $@

smrt100-2000.v2.best2.tsv:smrt100-2000.v2.tsv
	@../select-best.py 2 $< > $@

smrt100-2000-first.fna:smrt100-2000.fna
	@../../../scripts/select-seq.py -f 7 $< > $@

smrt100-2000-last.fna:smrt100-2000.fna
	@../../../scripts/select-seq.py -l 8 $< > $@

query-protein_f5.fsa:query-protein.fasta
	@../../../scripts/select-seq.py -f 5 $< > $@

Q6GZS4.fsa:ref-protein.fsa
	@../../../scripts/select-seq.py -p Q6GZS4 $< > $@

ref-query-protein-v1.tsv:ref-query-protein.tsv
	@cat $< | awk '/^[^#]/ {print $$1 "\t" $$2 "\t" $$3+$$4-1 "\t" $$5+$$6-1 "\t" $$7 "\t" $$8}' > $@
