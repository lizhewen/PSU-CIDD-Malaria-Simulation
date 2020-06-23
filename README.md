# PSU-CIDD-Malaria-Simulation

[Penn State](https://www.psu.edu/) - [Center for Infectious Disease Dynamics (CIDD)](https://www.huck.psu.edu/institutes-and-centers/center-for-infectious-disease-dynamics) - [Boni Lab](http://mol.ax/)

---

## Release Note

This version (v3.3) is developed for the Multiple-Drug-Resistant (MDR) paper. This codebase is based off Dang Nguyen Tran's v3.2 located in [this repo](https://github.com/maciekboni/PSU-CIDD-Malaria-Simulation).

This version has implemented Sporozoite Challenge codes, and a modified monthly reporter to output all mutation pair information throughout the simulation. **This version only supports single-location simulations.**

Under `mdr_analysis` folder, it contains several python scripts to parse output files, compute interquartile ranges for a series (default 100) of files from the same simulation settings, and create visualizations. Refer to the Jupyter Notebook `ipynb` files (named with leading numbers) to view available MDR plots.

---

## Output Files Description
One tab separated values (TSV) data file and one comma separated values (CSV) are generated if the default modified monthly reporter is selected for use. The tables listed below show you how these two files are organized, with group separators indicated by the sentinel value `-1111`.

**monthlydata_*n*.txt** - Summary data for the model generated at the end of each simulated month.

| Block | Column Number | Description | Dataframe Header |
| ---|--- | --- | ---|
| Summary | 1 |  Model time, number of days elapsed | time_elapsed |
| | 2 | Model time, calendar date as system time | system_time |
| | 3 - 5 | Model time, calendar date (Year, Month, Day) | year / month / day |
| | 6 | Seasonal factor | seasonal_factor |
| | 7 | Treatment coverage, probability to be treated (0 - 1) | treatment_coverage_0_1 |
| | 8 | Treatment coverage, probability to be treated (0 - 10) | treatment_coverage_0_10 |
| | 9 | Population size | population |
| | 10 | Group separator | sep |
| EIR and PfPR | 11 | EIR by location per year | eir |
| | 12 | Group separator | sep |
| | 13 | Blood slide prevalence, PfPR (age 2-10) | bsp_2_10 |
| | 14 | Blood slide prevalence, PfPR (age<5) | bsp_0_5 |
| | 15 | Blood slide prevalence, PfPR all | blood_slide_prevalence |
| | 16 | Group separator                                        | sep                      |
| Infections by Location        | 17               | Monthly Number of new infections                       | monthly_new_infection    |
|                               | 18               | Group separator                                        | sep                      |
| Treatments by Location        | 19               | Monthly Number of treatments                           | monthly_new_treatment |
|                               | 20               | Group separator                                        | sep                      |
| Clinical Episodes by Location | 21               | Monthly Number of Clinical Episodes                    | monthly_clinical_episode |
|                               | 22               | Group separator                                        | sep                      |
| Genotype Frequency | 23-150 | Each Genotype's Current Frequency | *{each genotype's encoding}* |
|  | ... | See genotype frequency discussion |  |

**mutpair_*n*.txt** - Information about all the mutations occurred throughout the simulation.

| Column Number(s) | Description | Dataframe Header |
| --- | --- | --- |
| 1 | Model time, number of days elapsed | time |
| 2 | Genotype Number muted from | from |
| 3 | Genotype Number muted to | to |

---

## Genotype Information

To view a complete list of all the genotypes (their encoding and number) and how effective they can be killed by a certain drug, use this [csv file](https://github.com/lizhewen/PSU-CIDD-Malaria-Simulation/blob/v3.3_MDR_Eric/mdr_analysis/ef2020.txt).

Genotypes are encoded into the application via the input file which contains a YAML entry for `genotype_info`. This structure is organized as `loci` with the assoicated `alleles` that may be mutated during model execution. The following table outlines loci and alleles that are included in the sample configuration file.

**Outline of loci and allel encoding**

| Locus | Allele | Short Name | Description |
| --- | --- | --- | --- |
| PfCRT | K76 | K | Lumefantrine resistant; Amodiaquine sensitive |
| | 76T | T |  Lumefantrine sensitive; Amodiaquine resistant |
| PfMDR1 | N86 Y184, Single-Copy |  NY-- | Partial-Lumefantrine and Partial-Amodiaquine resistant |
| | 86Y Y184, Single-Copy | YY-- | Amodiaquine resistant |
| | N86 184F, Single-Copy | NF-- | Lumefantrine resistant |
| | 86Y 184F, Single-Copy | YF-- | Partial-Amodiaquine and Partial-Lumefantrine resistant |
| | N86 Y184, Double-Copy | NYNY | Partial-Lumefantrine and Partial-Amodiaquine resistant, higher strength; Mefloquine selects strongly for multi-copy |
| | 86Y Y184, Double-Copy | YYYY | Amodiaquine resistant, higher strength; Mefloquine selects strongly for multi-copy |
| | N86 184F, Double-Copy | NFNF | Lumefantrine resistant, higher strength; Mefloquine selects strongly for multi-copy |
| | 86Y 184F, Double-Copy | YFYF | Partial-Amodiaquine and Partial-Lumefantrine resistant, higher strength; Mefloquine selects strongly for multi-copy |
| K13 Propeller | C580 | C | Artemisinin sensitive |
| | 580Y | Y | Artemisinin resistant |
| Plasmepsin 2-3 | Plasmepsin 2-3 Single-Copy | 1 | Piperaquine sensitive |
| | Plasmepsin 2-3 Double-Copy | 2 | Piperaquine resistant |
| Hypothetical locus for multiple use | naïve | x | Experimental use in the model |
| | mutant | X | Experimental use in the model |

The short name field requires addtional note since genotype results generated by the model are based upon the short names. For example:

> KNY--C1x

Indicates that the parasite has the K76 allele from the pfcrt locus (K), N86 Y184 one copy of pfmdr1 (NY--), C580 from K13 Propeller locus (C), Plasmepsin 2-3 one copy (1), and is a naïve copy of the hypotehtical locus (x). 

---

## About *mdr_analysis* Folder

---

## Contact

To contact us, please visit our [lab's website](http://mol.ax/).