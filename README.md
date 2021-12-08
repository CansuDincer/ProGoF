
# ProGoF

ProGoF is Python-based script collection that can be run to analyse target-disease association data from Open Targets Platform and to annotate them through [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), [PanelApp](https://panelapp.genomicsengland.co.uk/) and [Sequence Ontology](http://www.sequenceontology.org/) databases for interested list of genes.

### Requirements

* Python : 3.8
  + pandas : 1.1.3
  + argparse : 1.1
  + requests : 2.24.0
  + xml-python : 0.4.3

### Run

The scripts can be run with the following inputs:

*source* : This is the source for which the retrieval will be done. It can be ClinVar somatic and germline *(eva_somatic and eva, respectively)*, gene2phenotype *(gene2phenotype)*, PanelAPP *(genomics_england)*, and Open Targets Genetics *(ot_genetics_portal)*.

*version* : This is the version of the data retrieved from Open Targets Platform.

*score* : This is the minimum (and equal) value of the genetic association score threshold.

*datap* : This represents the input path. If not specified, the current directory will be used.

*resultp* : This represents the output path. If not specified, the current directory will be used.

*Note: If ot_genetics_portal will be run, then [variant_summary.txt](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz) file should be download into input folder.*

#### Example running

To construct association data frames for version 21.04 of Gene2Phenotype without any threshold for genetic association score:

```
python3 scripts/01_association_evidence.py -source gene2phenotype -datap input/ -resultp output/ -version 21.04 -score 0.0
```

After constructing all association data frames for interested sources:

To prioritize the targets:

```
python3 scripts/02_target_prioritisation.py -datap input/ -resultp output/ -version 21.04
```

### Output


Output files will be numbered:

* 02 - Annotated association data frames for each source
* 03 - GoF/LoF filtered association data frames for each source
* 04 - Drug mode of mechanisms integrated association data frames for each source
* 05 - Merged filtered association data frame
* 06 - Merged filtered association data frame of validation set
* 07 - Merged filtered target data frame (includes tractability)

Files including statistics will be created for representing total number of targets or associations. 


### Contact

It is a product of PhD rotation study of Cansu Dincer under the supervision of [Dr Gosia Trynka](https://www.sanger.ac.uk/person/trynka-gosia/), [Dr Ian Dunham](https://www.ebi.ac.uk/about/people/ian-dunham) and Dr Mohd Karim from Open Targets, Wellcome Sanger Institute and European Bioinformatics Institute. 

For any problems or feedback on ProGoF, you can contact [here](mailto:cd7@sanger.ac.uk).

### Terms and Conditions

Copyright (c) 2021 Genome Research Ltd.

Author: Cansu Dincer cd7@sanger.ac.uk

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation; either version 3 of the License, or (at your option) any later 
version. 

This program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
details. 

You should have received a copy of the GNU General Public License along with 
this program. If not, see <http://www.gnu.org/licenses/>. 
