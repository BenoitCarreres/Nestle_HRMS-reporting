# Nestle_HRMS-reporting
Report generator for targeted HRMS standard-addition-based analysis 

AUTHOR: Benoit Manuel Carreres, for Société des Produits Nestlé SA

## Installation

- Software was tested for R 4.0.2
- Each CMD files should be changed to point to a valid R-script.exe file.
- R dependencies needs to be installed. To do that, execute the script called "00_installDeps.cmd"

## Execution

- Drag-and-drop the input files onto the "01_generateReport.cmd" executable. Each input file should be properly formated. For more details, please see the scientific article linked bellow.
- Output file will be generated in the same folder as the input with the same name and the extra suffix "_report". If more than one file is given another suffix "_m" will be added.

## Special notes

- /!\ Some cells may be empty. If this is the case, requirements to generate that value are most likely not met. For example, missing RSQ value is probably because there are less than two valid injected values to construct the linear model.
- /!\ It is prefferable to input multiple files if they are part of the same analysis, because only one file will be the output.

## Metadata requirements

Metadata is an important part of the process. This information should be contained in the first part of the XLSX file: lines starting with "#".
Any metadata can be given in the XLSX file will be reported into the final report.

- Few metadata are required for the script to run:
    - "Concentration Unit"        |   Character value indicating the units of given concentrations.
    - "Carry over max"            |   Threshold value in percent. Maximum acceptable peak area value. This value is relative to the median value in QC for the first injection level.
    - "QC Solvent Repeatability"  |   Threshold value in percent. Maximum value for the Coefficient of Variations (CV, aka RDS) in the linear model.
    - "QC Solvent Linearity"      |   Threshold value in percent. Maximum deviation of any QC value used in the linear model.
    - "Linearity"                 |   Threshold for the linear model Coefficient of determination RSQ (R²).
    - "QC sample recovery"        |
    - "Screening Cut-off"         |   Giving a value for this one is not mandatory. If a value is given, the report will be generated as a Screening report type.

