This repository will help accomplish two things:

1. Feature Detection using TOPPAS/OPENMS
2. Create an exclusion list for the QExactive Mass Spectrometer


## Running on your own data

We can run this on Binder, click it below

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/lfnothias/IODA_MS/targeted_draft?filepath=IODA_notebooks_welcome.ipynb)


You will need to do two things:

1. Open IODA_TOPPAS_mztab_generation.ipynb to run the feature finding to create the mztab file
2. Open IODA_exclusion_from_mztab.ipyn to create the exclusion file to download


## Testing

![Unit Test](https://github.com/lfnothias/IODA_MS/workflows/Unit%20Test/badge.svg)

Unit tests are run using github actions. To run them manually:

```make test-unit```

The actual tests are in the ```test``` directory. 
