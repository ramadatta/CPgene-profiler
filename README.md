# CPgene-profiler
### We want to find the CPgenes from a bacterial assembly and plot a heatmap of CP genes from all the assemblies based on the categories of carbapenamases
### First step is to blast the Assembly against CP gene and find the list of the contigs with CP genes
### Second, the CP genes are assigned to the Assembly file and a table is created
### From the table a heatmap is drawn

[![Build Status](https://travis-ci.org/tseemann/snippy.svg?branch=master)]()
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)]()
![Don't judge me](https://img.shields.io/badge/Language-Perl_5-steelblue.svg)

# CPgene-profiler
A tool that profiles a list of Carbapenamase gene presence-absence and reports co-carriage within a assembly

## Author
[Prakki Sai Rama Sridatta](https://twitter.com/prakki_rama)

## Synopsis

CPgene-profiler checks for a list of CarbaPenamase (CP) genes from a list of genome assemblies
provided in fasta file format. It reports the profile of all the CP genes available in the genome assemblies 
in the format of simple heatmap. Apart from this, it also reports the presence of co-carriage of CP genes 
within an assembly. Other assembly statistics such as N50, N90, Assembly Size are calculated and plots of 
length distribution of CP gene contigs from the list of assemblies are reported.

## Quick Start
