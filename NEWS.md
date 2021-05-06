# iSTOP 0.2.0

* Added a `NEWS.md` file to track changes to the package.
* Fixed issue where `locate_codons()` failed to locate some `ANA`, `TNT`, `CNC`
  and `GNG` codons depending on sequence context. This change does not affect
  typical iSTOP usage. (#8) 
* Addressed `tidyverse` deprecation warnings for `data_frame` and
  `as_data_frame`.
