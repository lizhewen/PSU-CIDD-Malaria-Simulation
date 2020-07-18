# Multiple-Drug-Resistance Analysis & Visualization

This repo contains all the Python codes used to categorize, visualize, and analyze multiple-drug-resistance (MDR) evolution trends.

The Python files (end in `.py`) contain all the plotting functions (for visualizations) and their supporting functions. The Jupyter Notebook files (end in `.ipynb`) are showing the outcome measure calculation functions and evolution trend visualizations for the Multiple-Drug-Resistance (MDR) Analysis, focusing on comparing the Multiple-Firstline-Therapy (MFT) with other strategies.

This program is open-sourced under MIT License. There is much freedom with MIT License, but you must include the original license here and properly cite this as your source.

## Legend Coloring

For the genotype evolution trend plots, there are two options of calculating drug-resistance strength -

- Option 1: summarizes / categorizes by each genotype's drug resistance. This options uses `resistant_strength_calc` function and the resistance is formatted as `a-b`, where `a` being  how many drugs the genotype is resistant to, and `b` being how many mutation events occurred with the genotype. Refer to the function def for more info
- Option 2: summarizes by each genotype's drug efficacy

| Resistance Strength        | Drug Efficacy | Legend Color            |
| -------------------------- | ------------- | ----------------------- |
| 0-0                        | [90,∞)        | \#32CD32, green         |
| *N/A*                      | [80,90)       | \#3497FF, blue          |
| 1-1                        | *N/A*         | \#FAD996, light orange  |
| 1-2                        | [70,80)       | \#FFAB00, medium orange |
| 1-3                        | *N/A*         | \#C76400, dark orange   |
| 2-2                        | [60,70)       | \#F88379, coral red     |
| 2-3                        | [0,60)        | \#FF1A00, medium red    |
| 2-4, most dangerous double | *N/A*         | \#000000, black         |
| 3-5, most dangerous triple | *N/A*         | #800080, purple         |

