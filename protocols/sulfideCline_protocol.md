### Protocol for analyzing sulfide content using the Cline Reagent.

**Materials**

- N,N-dimethyl-p-phenylenediamine sulfate
- FeCl3
- 50% HCl
- 2.6% Zinc acetate
- 1% Zinc acetate
- Sodium Sulfide nonahydrate
- 15ml Falcon tubes.

**Preparing Cline Reagent (CR)**

*High Reagent (40-250uM)*

- 1.6g N,N-dimethyl-p-phenylenediamine sulfate
- 2.4g FeCl3
- 50ml HCl
- 50ml dH2O
- Wrap bottle in dark tape to keep light out.

*Low Reagent (5-50uM)*

- 0.4g N,N-dimethyl-p-phenylenediamine sulfate
- 0.6g FeCl3
- 50ml HCl
- 50ml dH2O
- Wrap bottle in dark tape to keep light out.

**Prep zinc acetate solutions**

*10% zinc acetate*

- 11.96 g zinc acetate dihydrate
- 100 ml sterile nanopore dH2O

*1% zinc acetate*

- 90 ml sterile nanopore dH2O
- 10 ml 10% zinc acetate


**Preparing Standards**

The bottle should be prepared under anaerobic conditions.

Stock solution: weight out ~100 mg of Na2S 9H2O, add to 10ml 1% Zinc acetate. Calculate the concentration of this, then use this to make a 1mM solutions (dilute into 1% ZnAc).

Using this stock solution, make standards of the following concentrations.

High Curve standards [Sulfide] (uM)

- H_1. 250
- H_2. 200
- H_3. 150
- H_4. 100
- H_5. 50
- H_6. 20

From the original stock solution, prepare a 100uM stock solution for preparing the low curve standards.

Low Curve standards [Sulfide] (uM)

- L_1. 50
- L_2. 30
- L_3. 20
- L_4. 10
- L_5. 5
- L_6. 2

**Protocol**

Divide the samples up into two groups based on their expected concentrations.
Any sample expected to be above 40uM should be analyzed using the high curve and any sample expected to be below 40uM should be analyzed using the low curve.

*Set up*

- Add 1ml of of each standard to 1.5ml Eppendorf tubes.
- Also prepare 4 blanks, using 1% ZnAc.
- Set up 12 samples for analysis, adding 1 ml of sample to the respective 1.5 ml Eppendorf tubes.
    - MAKE SURE to mix the sample well before pipetting, since the Zn precipitates out the sulfide. Use larger diameter pipette tips if possible.
    - Set up the fourth sample in triplicate. Try to make sure that this sample will be one that is within the curve, don't pick one that will be below detection.
    - For the 7th sample, set up a second tube with only 900 µl, then add 100 µl of the corresponding standard. This will be the matrix spike.
    - At the end, set up 1 ml of the fourth standard (continuing check standard, CCV) and a blank (continuing check blank, CCB).
- This can be set up a day ahead of time.


*Generate a standard curve and DL*

- Add 80ul of Cline reagent to each standard and blank.
- Mix by inverting and vortexing. Do one sample every 45 seconds for high curve.
- Incubate for 4 hr.
- Turn on the spectrophotometer at least 30 minutes before analysis.
- Zero the spectrophotometer at 667nm using nanopore water
- Measure absorbance (1 sample every 45 seconds):
  - For low curve, transfer 1ml standard to cuvettes and measure absorbance.
  - For high curve, transfer 250ul of standard to cuvettes and add 750ul of ZnAc and measure absorbance.

*Data analysis*

Calculate the blank value by averaging the absorbance for each of the blanks. Subtract the blank value from each of the standards.

To generate the standard curve, plot the concentrations against the corrected absorbance values. Generate a linear regression model. The R^2 value should be >= 0.99.

Calculate the detection limit using the following equation:

DL (mg/L) = (residual standard deviation) * 3 / (slope of regression)

*Re-analyzing samples*

For the high curve samples:

- If the sample absorbance is above the curve, dilute it and remeasure the absorbance.
- If sample readings are below the lowest standard, but it's close to the bottom of the curve, analyze the remaining 750 µl of sample, diluted with 250 µl of 1% Zn Ac.
- If sample readings are below the lowest standard by a great deal or if the previous method didn't get the samples within the curve, re-analyze them using the low curve.

For the low curve samples:

- If sulfide concentrations are above the curve, re-analyze them using the high curve.
- If sulfide concentrations are below the lowest standard or calculated to be below the DL, flag them as <DL.

All of these measurements will be recorded in my lab notebook and transferred to my computer.


*Data analysis*

All the data analysis for the sulfide data from this project will be done in R.
The R scripts found at `code/cleaning_scripts/clean_sulfide_data.R` contain the needed scripts.
This file first generates a metadata sheet for both the water column and the incubation sulfide samples.
It then outlines a function to process the data.

At the end of an analysis, the measurements are recorded in the `dataRaw/waterChemistry/sulfide` folder, in a file named `sulfide_yyyymmdd.xlsx`, which was copied from `SULFIDE_DATA_TEMPLATE.xlsx`.
I call the function on the raw data file.
If it passes the test, it gets read into the `data_for_review` folder.
I inspect the report file and the raw data.
If it all passes, I'll copy the results file into the `good_data` folder.

In my R script, I can then run the last part of the script, which reads all the good data in, attaches the metadata, and saves out the appropriate files.
It will also read out a csv of the samples that need to be run yet.

##### References:

Cline, Joel D. 1969. “Spectrophotometric Determination of Hydrogen Sulfide.” Limnology and Oceanography, 454–58. doi:10.4319/lo.1969.14.3.0454.

Reese, B. K., Finneran, D. W., Mills, H. J., Zhu, M. X. & Morse, J. W. Examination and Refinement of the Determination of Aqueous Hydrogen Sulfide by the Methylene Blue Method. Aquat. Geochemistry 17, 567–582 (2011).

http://parts.igem.org/Cline_Reaction_Assay_Protocol

https://corn.org/wp-content/uploads/2009/12/S-10.pdf
