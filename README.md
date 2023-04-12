# Genetic effects of tissue-specific enhancers in schizophrenia and hypertrophic cardiomyopathy
## Doctoral thesis by Emanuele Felice Osimo for Imperial College London
The link to the thesis will be included once publicly deposited.
Contact me at eosimo at ic.ac.uk

### Chapter Chapter 4 - EP-WAS development, and its internal validation, both on UK Biobank datasets.



This repository contains the code for work contained in Chapter 4. This is a Nextflow pipeline, in which I performed an EP-WAS (enhancer-based GWAS): in other words I measured associations with schizophrenia of SNPs falling within \textsc{Neural significant enhancers} (please see Chapter 2 for a rationale and methods of developing \textit{enhancer-based} partitions); to do so, I have used both \textit{dominant}, recessive, and additive inheritance models. 
The EP-WAS associations were calculated in UK Biobank, where they were initially internally validated

Run this pipeline with:

```
nextflow run main.nf -profile servername
```


