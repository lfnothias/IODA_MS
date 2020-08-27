# IODA MS

This repository will help accomplish two things:

1. Feature Detection using TOPPAS/OPENMS
2. Create an exclusion list for the QExactive Mass Spectrometer


## Running on your own data

We can run this on Binder, click it below

[![Binder](https://mybinder.org/badge_logo.svg)](https://notebooks.gesis.org/binder/v2/gh/lfnothias/IODA_MS/MS2Planner_merge_w_master_mingdocker?urlpath=lab/tree/IODA_MS2Planner.ipynb)

You will need to do two things:

1. Open IODA_TOPPAS_mztab_generation.ipynb to run the feature finding to create the mztab file
2. Open IODA_exclusion_from_mztab.ipyn to create the exclusion file to download


## Testing

![Unit Test](https://github.com/lfnothias/IODA_MS/workflows/Unit%20Test/badge.svg)

Unit tests are run using github actions. To run them manually:

```make test-unit```

The actual tests are in the ```test``` directory.

## Running the Docker in Linux

Run the following:

`make build-standalone`

`make run_standalone_notebook`

Access the notebook at `localhost:9000`
