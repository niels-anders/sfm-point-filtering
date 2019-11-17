# -*- coding: utf-8 -*-
"""
Multi-core processor of laz serialization
"""

import multiprocessing as mp
import subprocess
import os
import sys


def worker(pyfile, lazfile):
    """ Execute python file """
    try:
        job = subprocess.Popen(['python', pyfile, lazfile])
        job.wait()
    except OSError as error:
        print (error, pyfile, lazfile)
    return lazfile


def callback(lazfile):
    """ Callback: Write to screen when finished """
    print (lazfile, 'done')


def apply_async_with_callback(pyfile, lazfiles):
    """ Distribute jobs over CPU pool """
    cpu = mp.cpu_count() - 1

    print ('# cpus: %d' % cpu)
    pool = mp.Pool(processes=cpu)
    for lazfile in lazfiles:
        pool.apply_async(worker, (pyfile, lazfile), callback=callback)
    pool.close()
    pool.join()


if __name__ == '__main__':
    PYFILE = 'LazTools.py'
    PATH = sys.argv[1]

    LIST_OF_FILES = os.listdir(PATH)
    LAZFILES = []

    for LAZFILE in LIST_OF_FILES:
        if (LAZFILE[-3:] == 'las') or (LAZFILE[-3:] == 'laz'):
            LAZFILES.append(PATH + LAZFILE)
    apply_async_with_callback(PYFILE, LAZFILES)
