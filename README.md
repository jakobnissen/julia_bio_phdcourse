# Jakob Nissen's Julia course, October 2021
This is a part of a larger introductory Julia course for biologists.
The part in this repository covers some important packages of the BioJulia ecosystem.

This repo contains slides in PDF formats, as well as exercises/solution pairs.
There is a complimentary `data` directory which is too big to be part of this repo - it will be distributed alongside with this repo.
Send me an email if you are from the future and need a link to the data directory.
For now, the data is stored at:
https://drive.google.com/file/d/1gnFxXslTu8CTOzco08l0m562teXm5cAs/view?usp=sharing

## Content:
1. Introduction to BioJulia, BioSequences and kmers
2. BioJulia parsers and FASTX
3. Sequence alignment with BioAlignments

Older archived content:
* Population genetics with PopGen
* High throughput alignments with XAM

There is no introduction to Julia, the students are expected to already know the basics of the language.
The slides are superficial, providing only enough information to solve the exercises.
Hence, it is not a proper introduction to bioinformatics as a discipline, either.
The purpose of this part of the course is to:

1. Use BioJulia and bioinformatics as a source of exercises and examples to
   teach Julia
2. Introduce the BioJulia ecosystem

## How to do exercises
Most of the questions can be approached in multiple different ways.
Some solutions will be simpler, and/or much more efficient than others.
I recommend not worrying too much about efficiency and focus on getting the job done.

If you don't finish them, don't worry, it doesn't matter.
These exercises will not be graded.
They are simply meant as practise to get to use the software to solve various small tasks.

Most of the challenge of the exercises is to figure out what to actually do.
Therefore, you need to be able to search for approaches on your own. Here are a few tips:

* If you have an object `x`, try just typing it into the REPL to see what it looks like
* Use `typeof(x)` to get its type.
* Get help about a type `T` or a function `f` with `?T` and `?f`
* Use `dump(x)` to see the memory layout of `x`, and `dump(typeof(x))` to see the general memory layout of objects of the same type.
* To see what you can do with a type, try checking `methodswith(typeof(x), supertypes=true)`
* If you want to find functions about e.g. pairwise alignment, try running `apropos("pairwise alignment")`.
  This will search for all functions and types where "pairwise alignment" is part of the docstring.
* Try to find the documentation of the packages you are using and skim them

Solutions for the exercises are provided if you get stuck.
They are intentionally NOT optimized for speed.
PLEASE don't look at the solutions until after you've tried to solve the problem yourself.
Given the exercises are not graded, you will gain nothing from doing so.

All packages necessary to solve these exercises should be in the `Project.toml` file in each exercise directory.
Instantiate that to install the required packages.
