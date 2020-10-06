# ms-jaric-2020-fish-reaction-to-seiche
This repository hosts the data processing and analysis of the fish reaction to seiche events in Lake Milada, Czechia.

## Structure

```bash
├── ms-jaric-2020-fish-reaction-to-seiche.Rproj
|
├── data
│   ├── products    - any data derived from raw data
│   └── raw         - raw data
│
├── main.R          - main script
├── outputs         - any results
├── R               - scripts sourced by main.R
│
└── README.md
```

## How to start


1. Copy file `.Renviron_template` to a new file called `.Renviron` and update your database credentials there. The file is git ignored so the credentials should not be accidentaly pushed to the upstream repository. Finally, restart your R session - changes to .Renviron take effect when the new session is initialized.
1. Open `.Rproj` in Rstudio
1. Install package `renv`

   ```R
   install.packages("renv")
   ```

1. Install all required packages locally into new environment

   ```R
   renv::restore()
   ```

1. Run `main.R` script
