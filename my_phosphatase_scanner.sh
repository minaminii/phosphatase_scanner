#!/bin/bash

hmmscan -o GiardiaDB.general.hmmscan.out --tblout GiardiaDB.general.hmmscan.tab --domtblout Gi
ardiaDB.general.hmmscan.dom.tab ../phosphatase_domains.hmm3 GiardiaDB-3.1_GintestinalisAssemblageA_AnnotatedProteins.fasta

