from jug import TaskGenerator
import pandas as pd
from mirror import mirror_all_files, build_link_structure, create_ena_file_map
from links import build_links_sample
from config import MIRROR_BASEDIR, USE_ASPERA
import sys

build_links_sample = TaskGenerator(build_links_sample)

@TaskGenerator
def create_mirror(study_accession, target_directory):
    import ena
    filetable = ena.get_project_reads_table(study_accession, as_pandas_DataFrame=True)
    filetable = ena.expand_fastq_columns(filetable)
    mirror_all_files(filetable, MIRROR_BASEDIR, use_aspera=USE_ASPERA)
    if target_directory != '*':
        try:
            build_link_structure(filetable, MIRROR_BASEDIR, target_directory, study_accession)
        except IndexError:
            print("Failed to create link structure for:", study_accession, file=sys.stderr)
            raise
    return filetable

@TaskGenerator
def mirror_sample(sample_accession):
    import ena
    filetable = ena.get_project_reads_table(sample_accession, as_pandas_DataFrame=True)
    filetable = ena.expand_fastq_columns(filetable)
    mirror_all_files(filetable, MIRROR_BASEDIR, use_aspera=USE_ASPERA)
    return filetable

m = {}
studies = pd.read_table("studies.txt", comment='#')
for _, row in list(studies.iterrows()):
    target_directory = ('../' + row.directory_name if row.directory_name != '*' else '*')
    m[row.study_accession] = create_mirror(row.study_accession, target_directory)

create_ena_file_map(m, "filemap.tsv", MIRROR_BASEDIR)

links = []
samples = []
for s in open('samples.txt'):
    samples.append(mirror_sample(s.rstrip()))
    links.append(build_links_sample(samples[-1]))
