from os import symlink, makedirs, path

from config import MIRROR_BASEDIR, TARGET_BASEDIR

def build_links_sample(filetable):
    if len(set(filetable.sample_accession.values)) > 1:
        raise ValueError("Expected just 1 sample???")
    [s] = set(filetable.sample_accession.values)

    files = filetable.ftp.map(lambda n : n[len('ftp.sra.ebi.ac.uk/'):]).values

    makedirs(f'{TARGET_BASEDIR}/{s}', exist_ok=True)
    for f in files:
        symlink(f'{MIRROR_BASEDIR}/{f}',f'{TARGET_BASEDIR}/{s}/'+path.basename(f))
