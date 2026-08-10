# Bugs fixed relative to upstream RIsearch 1.2

Defects in upstream [RTH-tools/risearch](https://github.com/RTH-tools/risearch) that
this fork fixes. Each was reproduced against the original C and has a regression test.

* Fixed infinite loop with FASTA input files where the last line is a header
* A record with no sequence, or one of only gap characters, will overflow
* Alignment with no complementary pair would access invalid memory
* Corrupt bytes in a sequence were dropped without a word, so a damaged record was
  read as a shorter one that still looked valid, and scored as if it were whole
* Fixed isalpha (undefined behavior for bytes above 127, since a plain char is signed)
* Switched isalpha with a table lookup, improved correctness and FASTA reading speeds
* Leaked memory when sequence was rejected `-Q "ACGU@ACGU"`

Behavioral changes:

* Sequence with a control byte, or a byte of value 127 or above, is rejected instead
  of omitted
