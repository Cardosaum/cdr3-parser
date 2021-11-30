# cdr3-parser

A simple program to convert `aafreq` files (generated as an intermediate step
from [ATTILA](https://github.com/Cardosaum/attila)) into `csv` or `json` files
containing statistics about all VH-CDR3 regions presented in the input file.

## Overview

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

## Installation

The project is still in 'alpha', only intended for local use at the University
of Brasília's Bioinformatics Laboratory, but if you want to run it locally, you
can compile the source code. Just bear in mind that you'll need a working rust
installation in order to do that.

    git clone https://github.com/Cardosaum/cdr3-parser.git
    cd cdr3-parser
    cargo install --path .

After successfully compiling the code, and assuring that `~/.cargo/bin` is in
your `$PATH`, you are good to go :)
    
