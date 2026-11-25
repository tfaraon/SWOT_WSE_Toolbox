# Contributing

Open an issue with the product version, raster spacing, station settings and the smallest example that reproduces the problem. A synthetic granule is sufficient. Remove Earthdata credentials from logs and notebooks before sharing them.

For a change to a filter, describe which observations it affects. Include a test with a known expected value and run `python -m unittest discover -s tests -v` from the repository root. Keep discussions focused on the analysis and treat contributors respectfully.
