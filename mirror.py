import os
from os import path, makedirs
import pandas as pd
import urllib
import pathlib
import requests
from contextlib import closing
from math import isnan
from jug import bvalue
from config import ASPERA_BINARY, ASPERA_KEY


def mirror_path(mirror_basedir, ftp):
    import urllib
    import pathlib
    url = urllib.parse.urlparse('http://' + ftp)
    p = pathlib.PurePath(url.path)
    target_dir = mirror_basedir / p.parent.relative_to('/')
    return target_dir / p.name



def md5sum_file(ifile):
    '''Computes MD5 sum of ifile'''
    import hashlib
    m = hashlib.md5()
    BLOCK_SIZE = 8192
    with open(ifile, 'rb') as ifile:
        while True:
            data = ifile.read(BLOCK_SIZE)
            if not data:
                return m.hexdigest()
            m.update(data)

def http_download_file(url, ofile):
    with closing(requests.get(url, stream=True)) as ifile, \
                open(ofile, 'wb') as ofile:
        for chunk in ifile.iter_content(8192):
            ofile.write(chunk)

def aspera_download_file(aspera_url, ofile):
    '''Call ascp on the command line to download `aspera_url` to `ofile`'''
    import subprocess
    cmdline = [
            ASPERA_BINARY,
            '-P33001', # Use special port
            '-T', # No encryption
            '-l', '300m',
            '-i', ASPERA_KEY,
            aspera_url,
            str(ofile)]
    subprocess.run(cmdline, check=True)

def mirror_all_files(filetable, mirror_basedir, *, progress=True, use_aspera=False):
    n = len(filetable)
    for i in range(n):
        if progress:
            print("Processing file {} of {}.".format(i + 1, n))
        source = filetable.iloc[i]
        if type(source.ftp) == float and isnan(source.ftp):
            continue

        urlraw = 'http://' + source.ftp
        url = urllib.parse.urlparse(urlraw)
        p = pathlib.PurePath(url.path)
        target_dir = mirror_basedir / p.parent.relative_to('/')

        makedirs(target_dir, exist_ok=True)
        ofile = target_dir / p.name

        if path.exists(ofile):
            if os.stat(ofile).st_size != int(source.bytes):
                print("Existing output file has wrong size. Removing...")
                os.unlink(ofile)
#            elif md5sum_file(ofile) != source.md5:
#                print("Existing output file has wrong hash. Removing...")
#                os.unlink(ofile)
            else:
                print("Correct output file exists. Skipping...")
                continue
        if use_aspera:
            aspera_url = 'era-fasp@fasp.sra.ebi.ac.uk:'+url.path
            aspera_download_file(aspera_url, ofile)
        else:
            http_download_file(urlraw, ofile)


def norm_path(p):
    p = str(p)
    if p.endswith('_1.fastq.gz'):
        return pathlib.PurePath(p[:-len('_1.fastq.gz')]+'.pair.1.fq.gz')
    if p.endswith('_2.fastq.gz'):
        return pathlib.PurePath(p[:-len('_2.fastq.gz')]+'.pair.2.fq.gz')
    if p.endswith('.fastq.gz'):
        return pathlib.PurePath(p[:-len('.fastq.gz')] + '.single.fq.gz')
    raise ValueError("Cannot normalize {}".format(p))


def build_link_structure(filetable, mirror_basedir, data_basedir, sample_fname):
    data_basedir = pathlib.PurePath(data_basedir)
    makedirs(data_basedir, exist_ok=True)

    with open(data_basedir / sample_fname, 'w') as samplefile:
        for s in set(filetable.sample_accession):
            samplefile.write("{}\n".format(s))

    prefix_fields = [ col for col in 
        ('library_layout', 'library_strategy', 'library_source', 'library_selection')
        if filetable[col].nunique() > 1
    ] 

    n = len(filetable)
    for i in range(n):
        source = filetable.iloc[i]
        if type(source.ftp) == float and isnan(source.ftp):
            continue

        target = data_basedir

        if prefix_fields:
            target = target / "_".join( source[s] for s in prefix_fields ) 

        target = target / source.sample_accession 
        makedirs(target, exist_ok=True)
        target = target / norm_path(source.ftp).name

        os.symlink(mirror_path(mirror_basedir, source.ftp), target)

def create_ena_file_map(studies_tables, vol_map, MIRROR_BASEDIR):
    def annotate_link(p):
        def drop_hostname(addr):
            while not addr.startswith("vol1"):
                newaddr = addr.split('/', 1)
                if len(newaddr) == 1:
                    raise ValueError("Couldn't find mirror root")
                addr = newaddr[-1]
            return addr

        if pd.isnull(p):
            return (p, p)

        p = drop_hostname(p)
        if p.endswith('_1.fastq.gz'):
            return ("fastq_1", p)
        if p.endswith('_2.fastq.gz'):
            return ("fastq_2", p)
        if p.endswith('.fastq.gz'):
            return ("fastq_single", p)
        raise ValueError("Cannot annotate {}".format(p))

    with open(path.join(MIRROR_BASEDIR, vol_map), 'w') as out:
        out.write("#study_accession\trun_accession\tsample_accession\texperiment_accession\tfastq_1\tfastq_2\tfastq_single\n")
        for study in studies_tables:
            data = bvalue(studies_tables[study])  ## bvalue because this is a Tasklet
            if data is None:
                print("We got no file information for", study, ". Incomplete/older jug internal state? Skipping", study)
                continue
            annotations = data["ftp"].map(annotate_link)
            data["filetype"], data["filepath"] = zip(*annotations)

            subdata = data[["study_accession", "run_accession", "sample_accession",
                            "experiment_accession"]].drop_duplicates()
            files = data.pivot(values="filepath", columns="filetype")
            grouped = pd.merge(subdata, files, left_index=True, right_index=True)

            for _, record in grouped.iterrows():
                # 6 columns: project_accession, sample_accession, experiment_accession, fastq_1, fastq_2, fastq_single
                try:
                    fastq_1 = record.fastq_1
                except:
                    fastq_1 = ''
                else:
                    if fastq_1 is None:
                        fastq_1 = ''
                    else:
                        fastq_1 = path.join(MIRROR_BASEDIR, fastq_1)

                try:
                    fastq_2 = record.fastq_2
                except:
                    fastq_2 = ''
                else:
                    if fastq_2 is None:
                        fastq_2 = ''
                    else:
                        fastq_2 = path.join(MIRROR_BASEDIR, fastq_2)

                try:
                    fastq_single = record.fastq_single
                except:
                    fastq_single = ''
                else:
                    if fastq_single is None:
                        fastq_single = ''
                    else:
                        fastq_single = path.join(MIRROR_BASEDIR, fastq_single)

                out.write(f"{record.study_accession}\t{record.run_accession}\t{record.sample_accession}\t{record.experiment_accession}\t{fastq_1}\t{fastq_2}\t{fastq_single}\n")

