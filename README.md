# Mirror selected ENA data locally

This was developed to support the following publication:

> Coelho, L.P., Alves, R., del Río, Á.R. et al. Towards the biogeography of
> prokaryotic genes. Nature 601, 252–256 (2022).
> https://doi.org/10.1038/s41586-021-04233-4

Please cite this publication if you use this tool.

To use, create a TSV file called `studies.txt` with the following columns
(first row should be column headers):

    study_accession	directory_name	reference	comment

Then, run `jug execute` to retrieve all the data. If you later add more rows,
Jug will ensure that only missing data will be downloaded.

The file `config.py` adds a fews configuration options. You can start with the
configuration in `config-example.py` and edit it to suit your needs. In
particular, using aspera, which is turned off by default, can significantly
speed up downloads.

Dependencies

- pandas
- [requests](http://docs.python-requests.org/en/master/)
- [jug](http://jug.rtfd.io)

