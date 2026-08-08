# RIsearch1 (tauso fork)

RIsearch1: RNA–RNA, RNA-DNA and DNA-DNA interaction prediction using a simplified
nearest-neighbor energy model.

This fork carries **RIsearch1 only**. RIsearch2 and the siRNA off-target pipeline
live upstream at [RTH-tools/risearch](https://github.com/RTH-tools/risearch); they
were removed here because nothing downstream of this fork uses them. Their history
is still in this repository if you need it.

## Python package

A precompiled binary from this fork is published to PyPI as `risearch-tauso`, so
RIsearch is a normal Python dependency for tauso and any other downstream tool:

```bash
pip install risearch-tauso
```

```python
import risearch_tauso, subprocess
subprocess.run([risearch_tauso.executable_path(), "-q", "query.fa", "-t", "target.fa"])
```

Or as a CLI shim, which forwards all arguments straight to the bundled binary:

```bash
risearch-tauso -q query.fa -t target.fa
python -m risearch_tauso -q query.fa -t target.fa
```

The PyPI package is **not** canonical upstream RIsearch — it is the tauso-team
fork. Use upstream if you want the unmodified tool.

## Building from source

Requires a C++17 compiler and CMake 3.20 or newer.

```bash
cmake -S RIsearch1 -B RIsearch1/build -DCMAKE_BUILD_TYPE=Release
cmake --build RIsearch1/build -j
```

That produces `RIsearch1/bin/RIsearch` and `RIsearch1/bin/RIsearch.dbg`, the
second with the debug tracing compiled in. The wheel build runs the same CMake via
`setup.py`.

## Tests

```bash
./RIsearch1/build/tests/risearch1_tests
```

Unit tests for the nucleotide coding, the min/max helpers, the alignment symbols,
the energy matrix and the FASTA reader; end-to-end tests that run `main()` in
process with a constructed argv and assert on its output. Pass
`-DRISEARCH_BUILD_TESTS=OFF` to skip googletest entirely, as the wheel build does.

## Running

```bash
RIsearch -q query.fa -t target.fa
```

Both files may hold several sequences; RIsearch scans all against all. Single
sequences can be given directly with `-Q acgu -T acgu`. See `RIsearch1/Manual.pdf`
for the full option list.

## Copyright

Copyright 2021 by the contributors; see `RIsearch1/README`.

RIsearch1 is released under the GNU General Public License version 3. This is free
software: you can redistribute it and/or modify it under the terms of that licence,
either version 3 or (at your option) any later version. You should have received a
copy of the GNU General Public License along with RIsearch — see the file COPYING.
If not, see <http://www.gnu.org/licenses/>.

This software is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Citation

If you use RIsearch in a publication, please cite:

**RIsearch: fast RNA-RNA interaction search using a simplified nearest-neighbor
energy model.** Wenzel A, Akbasli E, Gorodkin J. *Bioinformatics*. 2012 Nov
1;28(21):2738-46.

## Contact

For problems with this fork, open an issue here. For upstream RIsearch:
<software+crispron@rth.dk>
