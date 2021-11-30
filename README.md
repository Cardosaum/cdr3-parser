# cdr3-parser

A simple program to convert `aafreq` files (generated as an intermediate step
from [ATTILA](https://github.com/Cardosaum/attila)) into `csv` or `json` files
containing statistics about all VH-CDR3 regions presented in the input file.

```bash
$ cdr3-parser --help

cdr3-parser 0.1.0
Output chemical properties and statistics of CDR3 sequences within a `aafreq` file.
You may choose either `json` or `csv` as output formats.

USAGE:
    cdr3-parser [FLAGS] <INPUT_FILE>

FLAGS:
        --cdr-only    Informs program that the input contains only CDR3VH sequences, separated by a new line
    -c, --csv         Select `csv` as output format [default]
    -h, --help        Prints help information
    -j, --json        Select `json` as output format
    -V, --version     Prints version information

ARGS:
    <INPUT_FILE>    Sets the input file to use
```
