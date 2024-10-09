# GTTL Based Implementation of the NTCard Algorithm

This tool is an implementation of the
[NTCard](https://doi.org/10.1093%2Fbioinformatics%2Fbtw832)
algorithm by Mohamadi, Khan and Birol, which estimates the k-mer count in a
given sequence file.

## Advantages

This implementation has the following advantages over the original implementation
given in the aforementioned paper:

- `ntcard_mn.x` is 7-8x faster (on short-read data)
- `ntcard_mn.x` is generally more accurate than the original `ntcard`.
    - On larger short-read samples, the relative error generally remains <1%.
- `ntcard_mn.x` can read `gzip`-compressed files directly.
- `ntcard_mn.x` can use more than one thread if necessary.
    - This comes at the cost of virtual memory, since the file will be `mmap`ed
    - This is only ever effective on files that are not `gzip` compressed, since decoding 
    from the `gzip` file will always be single-threaded.

## Building

Calling `make` in this directory will build the tool, assuming that the
depedency on GTTL and its dependencies are satisfied and the variable `GTTL`
in the `Makefile` is set accordingly.
  
An executable binary is then generated under the name `ntcard_mn.x`.

## Usage

- The `-b` flag may significantly reduce runtime and memory-usage.
- `ntcard_mn.x` can handle files of the `fastq` and `fasta` formats, as well as their
  `gzip`-compressed forms.


For further usage information, try `./ntcard_mn.x --help`

### Examples

For instance, the F0-values (ie. the number of distinct k-mers in the sample)
for `k=20` can be determined as follows:
`./ntcard_mn.x -q 20 -b <inputfile.fastq.gz>`
