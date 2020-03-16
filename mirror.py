import os
import shutil
import urllib
import pathlib
import requests
import pandas as pd
from jug import bvalue
from math import isnan
from contextlib import closing
from os import path, makedirs
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


def wget_download_file(url, ofile):
    '''Call wget on the command line to download `url` to `ofile`'''
    import subprocess
    cmdline = ['wget',
               url,
               '-O',
               str(ofile)]
    print('WGET_CMD', cmdline)
    subprocess.run(cmdline, check=True)


def aspera_download_file_temp_dir(aspera_url, ofile):
    '''Call ascp on the command line to download `aspera_url` to `ofile`
     - Temporaily downloads to bork9 first and then copies to ofile due to
    data fragmentation issies encountered on scb2.'''
    import subprocess
    temp_download_area = '/g/bork9/fullam/temp_downloads/'
    temp_download_name = os.path.join(temp_download_area,
                                      os.path.basename(ofile))
    cmdline = [ASPERA_BINARY,
               '-P33001',  # Use special port
               '-T',  # No encryption
               '-l', '300m',
               '-i', ASPERA_KEY,
               aspera_url,
               str(temp_download_name)]
    print('ASPERA_CMD', cmdline)
    print('DESTINATION: ', ofile)
    subprocess.run(cmdline, check=True)
    if os.path.isfile(temp_download_name):
        shutil.move(temp_download_name, ofile)


def aspera_download_file(aspera_url, ofile):
    '''Call ascp on the command line to download `aspera_url` to `ofile`'''
    import subprocess
    cmdline = [ASPERA_BINARY,
               '-P33001',  # Use special port
               '-T',  # No encryption
               '-l', '300m',
               '-i', ASPERA_KEY,
               aspera_url,
               str(ofile)]
    print('ASPERA_CMD', cmdline)
    subprocess.run(cmdline, check=True)


def mirror_all_files(study_accession, filetable, mirror_basedir, *, progress=True, use='HTTP'):
    n = len(filetable)
    for i in range(n):
        if progress:
            print("Processing file {} of {}. Study: {}".format(i + 1, n, study_accession))
        source = filetable.iloc[i]
        if type(source.ftp) == float and isnan(source.ftp):
            continue

        urlraw = 'http://' + source.ftp
        url = urllib.parse.urlparse(urlraw, allow_fragments=False)
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
        if use == 'ASPERA':
            aspera_url = 'era-fasp@fasp.sra.ebi.ac.uk:' + url.path
            for attempt_number in range(3):
                aspera_download_file(aspera_url, ofile)
                if check_file(ofile, source):
                    break
        elif use == 'WGET':
            for attempt_number in range(3):
                wget_download_file(source.ftp, ofile)
                if check_file(ofile, source):
                    break
        else:
	    for attempt_number in range(3):
                http_download_file(urlraw, ofile)
                if check_file(ofile, source):
                    break


def check_file(ofile, source):
    print('Checking file...')
    file_size = os.stat(ofile).st_size
    print('File size: {0}'.format(file_size))
    if file_size != int(source.bytes):
        print("File has wrong size. {0} vs.{1}".format(file_size,
                                                       source.bytes))
        print('Removing..')
        os.unlink(ofile)
        return False
    file_md5sum = md5sum_file(ofile)
    print('File md5: {0}'.format(file_md5sum))
    if file_md5sum != source.md5:
        print("File has wrong md5. {0} vs.{1}".format(file_md5sum,
                                                      source.md5))
        print('Removing..')
        os.unlink(ofile)
        return False
    print("Downloaded file OK...")
    return True


def norm_path(p):
    p = str(p)
    if p.endswith('_1.fastq.gz'):
        return pathlib.PurePath(p[:-len('_1.fastq.gz')] + '.pair.1.fq.gz')
    if p.endswith('_2.fastq.gz'):
        return pathlib.PurePath(p[:-len('_2.fastq.gz')] + '.pair.2.fq.gz')
    if p.endswith('.fastq.gz'):
        return pathlib.PurePath(p[:-len('.fastq.gz')] + '.single.fq.gz')
    # if p.endswith('.fq1.gz'):
    #     return pathlib.PurePath(p[:-len('.fq1.gz')] + '.single.fq.gz')
    raise ValueError("Cannot normalize {}".format(p))


def build_link_structure(filetable, mirror_basedir, data_basedir, sample_fname):
    data_basedir = pathlib.PurePath(data_basedir)
    makedirs(data_basedir, exist_ok=True)

    with open(data_basedir / sample_fname, 'w') as samplefile:
        for s in set(filetable.sample_accession):
            samplefile.write("{}\n".format(s))

    prefix_fields = [col for col in
                     ('library_layout',
                      'library_strategy',
                      'library_source',
                      'library_selection')
                     if filetable[col].nunique() > 1]

    n = len(filetable)
    for i in range(n):
        source = filetable.iloc[i]
        if type(source.ftp) == float and isnan(source.ftp):
            continue

        target = data_basedir

        if prefix_fields:
            target = target / "_".join(source[s] for s in prefix_fields)

        target = target / source.sample_accession
        makedirs(target, exist_ok=True)
        target = target / norm_path(source.ftp).name

        try:
            os.symlink(mirror_path(mirror_basedir, source.ftp), target)
        except FileExistsError as e:
            os.unlink(target)
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
        if p.endswith('R1.fastq.gz'):
            return ("fastq_1", p)
        if p.endswith('_2.fastq.gz'):
            return ("fastq_2", p)
        if p.endswith('R2.fastq.gz'):
            return ("fastq_2", p)
        if p.endswith('.fastq.gz'):
            return ("fastq_single", p)
        # else:
        #     return ("fastq_single", p)
        raise ValueError("Cannot annotate {}".format(p))

    with open(path.join(MIRROR_BASEDIR, vol_map), 'w') as out:
        out.write("#study_accession\trun_accession\tsample_accession\texperiment_accession\tfastq_1\tfastq_2\tfastq_single\n")
        for study in studies_tables:
            data = bvalue(studies_tables[study])  ## bvalue because this is a Tasklet
            if data is None or len(data) == 0:
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
