# Mirror selected ENA data locally

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

