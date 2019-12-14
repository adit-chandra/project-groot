# Project Groot
 <img src="etc/groot.png" width="800">

*An Algorithmic Approach to Phylogenic Inference on Reduced Genomic Data*

### Dependecies
- `python 3.7`
- binaries are compiled for `Ubuntu 16.04.6 LTS`

### Usage
To run inference procedure, execute:

`python3 groot.py -d <data_dir> -mem <mem_size> -t <n_threads> -k <klen> --ks <kslen> -n <filter> -l <long>`

See `groot.py` for full option descriptions.

**NOTE**: `data_dir` must contain sequences in [FASTA](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp) format.

### References
- Baum, D. A., Smith, S. D. (2013). Tree Thinking: An Introduction to Phylogenetic Biology. Greenwood Village, CO: Roberts and Company.
- Felsenstein J (2004). Inferring Phylogenies. Sunderland, Massachusetts: Sinauer Associates. ISBN 978-0-87893-177-4.
- Sankoff D, Morel C, Cedergren RJ (October 1973). ”Evolution of 5S RNA and the non-randomness of base replacement”. Nature. 245 (147): 232–4. doi:10.1038/newbio245232a0. PMID 4201431
- Wheeler WC, Gladstein DS (1994). ”MALIGN: a multiple nucleic acid sequence alignment program”. Journal of Heredity. 85 (5): 417–418. doi:10.1093/oxfordjournals.jhered.a111492.
- Cruaud, A. et al (2014). Empirical assessment of RAD sequencing for interspecific phylogeny. Mol. Biol. Evol. 31, 1272–1274.
- Andrews KR, Good JM, Miller MR, Luikart G, Hohenlohe PA (2016) Harnessing the power of RADseq for ecological and evolutionary genomics. Nature Reviews Genetics, 17, 81–92.
- Chong Z, Ruan J, Wu C-I (2012) Rainbow: an integrated tool for efficient clustering and assembling RAD-seq reads. Bioinformatics, 28, 2732–2737.
- Fan H, Ives AR, Surget-Groba Y, Cannon CH (2015) An assembly and alignment-free method of phylogeny reconstruction from next-generation sequencing data. BMC Genomics, 16(1), 349.
- Stuart GW, Moffett K, Leader JJ. A Comprehensive Vertebrate Phylogeny Using Vector Represen- tations of Protein Sequences from Whole Genomes. Mol Biol Evol. 2002;19:554–62.
- Prasad AB, Allard MW, NISC Comparative Sequencing Program, Green ED. Confirming the Phylogeny of Mammals by Use of Large Comparative Sequence Data Sets. Mol Biol Evol. 2008;25:1795–808.
