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

1. Open `.Rproj` in Rstudio
2. Install package `renv`

   ```R
   install.packages("renv")
   ```

3. Install all required packages locally into new environment

   ```R
   renv::restore()
   ```

4. Copy file `.Renviron_template` into `.Renviron` file and update your database credentials there. The file is git ignored so the credentials should not be accidentaly pushed to public repo.

5. Run `main.R` script
