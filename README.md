# KruvelabFns
Functions that Kruvelab uses on regular basis.


To install the R-package run:
```
devtools::install_github("kruvelab/kruvelab_functions")
```

## Functions overview

All functions developed or currently in progress can be found [here](https://kruvelab-my.sharepoint.com/:x:/g/personal/idarahusu_kruvelab_onmicrosoft_com/ESnA3VMZ9lBNubtbL7UYJlYBBiDkiMGXWNy8Pc06lljlNA?e=fkE42b).
In addition to the status of the function, one can find information about inputs, outputs, dependencies, and the person in charge of the function.

## Developed Google colab workbooks for functions in Python

Here you can find workbooks developed for applying specific functions in Python.

[Workflow for processing data in SMILES format](https://colab.research.google.com/drive/1PM5D8NqRcGR2x5KlWvjKf1gyky1OCv1H?usp=sharing)
* PART A. File formats
* PART B. Standardization of SMILES notations
* PART C. Calculation of Mordred descriptors
* PART D. Visualization of chemical space


## When you are developing a function...

* use lowercase letters and underscores for function names and variables (except for names such as SIRIUS, etc., or using one variable to calculate other(s) such as "SMILES2descriptors").
* check if the same variables are used in other functions; if yes, then keep to the same writing of the variable.
* specify variable types for inputs if possible.
* note down which libraries does your function depend on.
* use a [function_template](https://github.com/kruvelab/kruvelab_functions/blob/main/in_work/function_template.R) to document your function.

