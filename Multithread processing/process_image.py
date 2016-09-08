#! /usr/bin/env python
#
# Requirements:
#
# Python 3
# pip3 install pyfits
# pip3 install astropy
# pip3 install astroml
# pip3 install sklearn

import warnings

warnings.filterwarnings("ignore")

import sys, os, getopt
import os.path
import time
import numpy as np
import pyfits
import multiprocessing as mp
import queue
import ctypes

from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from astroML.plotting import setup_text_plots

setup_text_plots(fontsize=8, usetex=True)

sigmatofwhm = 2 * np.sqrt(2 * np.log(2))
fwhmtosigma = 1. / sigmatofwhm

class Worker_convolve(mp.Process):
    def __init__(self, work_queue, result_queue):

        # base class initialization
        mp.Process.__init__(self)

        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False

    def run(self):
        while not self.kill_received:

            # get a task
            try:
                x_range, y_range = self.work_queue.get_nowait()
            except queue.Empty:
                break

            # the actual processing
            print('Running process id: ', mp.current_process().name)

            im_size = np.asarray(shared_im.shape)
            kernel_size = np.asarray(shared_kernel.shape)

            im_pad = (kernel_size - 1) / 2 + 10
            x_pad = [0 if x_range[0] == 0 else -im_pad[0], 0 if x_range[1] == (im_size[0] - 1) else im_pad[0]]
            y_pad = [0 if y_range[0] == 0 else -im_pad[1], 0 if y_range[1] == (im_size[1] - 1) else im_pad[1]]

            shared_nim[x_range[0]:x_range[1], y_range[0]:y_range[1]] = convolve(
                shared_im[x_range[0] + x_pad[0]:x_range[1] + x_pad[1], y_range[0] + y_pad[0]:y_range[1] + y_pad[1]],
                shared_kernel, normalize_kernel=True)[-x_pad[0]:x_range[1] - x_range[0] - x_pad[0],
                                                                       -y_pad[0]:y_range[1] - y_range[0] - y_pad[0]]

            #print('Worker_convolve done: ', mp.current_process())
            self.result_queue.put(id)

class Worker_median(mp.Process):
    def __init__(self, work_queue, result_queue):

        # base class initialization
        mp.Process.__init__(self)

        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False

    def run(self):
        while not self.kill_received:

            # get a task
            try:
                x_range, y_range = self.work_queue.get_nowait()
            except queue.Empty:
                break

            # the actual processing
            print('Running process id: ', mp.current_process().name)

            im_size = np.asarray(shared_im.shape)

            shared_nim[x_range[0]:x_range[1], y_range[0]:y_range[1]] = \
                shared_im[x_range[0]:x_range[1], y_range[0]:y_range[1]] - \
                np.nanmedian(shared_im[x_range[0]:x_range[1], y_range[0]:y_range[1]])

            #print('Worker_median done - id: ', mp.current_process().name)
            self.result_queue.put(id)


if __name__ == "__main__":

    argv=sys.argv[1:]
    im_file_in='data/VCC1010.MASKED.fits'
    do_recipe='median'
    do_overwrite=False
    n_cpu = 1
    n_core = 4
    # You can also ask the os how many processors has
    # n_cpu = mp.cpu_count()
    # n_core = 4 #Should be your physical threads perhaps?
    n_x=1
    n_y=6
    hyper_thread=1 # If you have hyper threading, then set the value to 2, otherwise set to 1

    try:
        opts, args = getopt.getopt(argv, "hfri:o:", ["recipe=", "ifile=", "ofile=","ncpu=","ncore=","nx=","ny="])
    except getopt.GetoptError:
        print('process_image.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('convolve_image.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            im_file_in = arg
        elif opt in ("-o", "--ofile"):
            im_file_out = arg
        elif opt in ("-r", "--recipe"):
            do_recipe = arg
        elif opt in ("--ncpu"):
            n_cpu = int(arg)
        elif opt in ("--ncore"):
            n_core = int(arg)
        elif opt in ("--nx"):
            n_x = int(arg)
        elif opt in ("--ny"):
            n_y = int(arg)
        elif opt == '-f':
            do_overwrite=True

    n_process = np.min( [n_cpu * n_core * hyper_thread, n_x*n_y] )
    grid_n = [n_y, n_x]  # Number of bins to divide image along [Y,X]
    im_file_out = 'data/VCC1010.MASKED.'+do_recipe.upper()+'.fits'

    print('Reading image: ', im_file_in)
    hdulist = pyfits.open(im_file_in)
    im_data = hdulist[0].data
    im_h = hdulist[0].header
    hdulist.close()
    im_size = np.asarray(im_data.shape)
    print('Image size: ', im_size)

    shared_im_base = mp.Array(ctypes.c_float, im_data.size)
    shared_im = np.ctypeslib.as_array(shared_im_base.get_obj())
    shared_im = shared_im.reshape(im_size[0], im_size[1])
    shared_im[:] = im_data

    shared_nim_base = mp.Array(ctypes.c_float, im_data.size)
    shared_nim = np.ctypeslib.as_array(shared_nim_base.get_obj())
    shared_nim = shared_nim.reshape(im_size[0], im_size[1])
    shared_nim[:] = 0.

    if os.path.isfile(im_file_out) is False or do_overwrite is True:

        tic = time.time()
        if do_recipe=='median':

            work_queue = mp.Queue()
            grid_n = np.asarray(grid_n)
            grid_mesh = np.ceil(im_size * 1. / grid_n).astype(int)
            grid_x = np.append(np.arange(0, im_size[0], grid_mesh[0]), im_size[0])
            grid_y = np.append(np.arange(0, im_size[1], grid_mesh[1]), im_size[1])

            for i in range(0, grid_x.size - 1):
                for j in range(0, grid_y.size - 1):
                    if work_queue.full():
                        print("Oh no! Queue is full after only %d iterations" % j)

                    x_range = [grid_x[i], grid_x[i + 1]]
                    y_range = [grid_y[j], grid_y[j + 1]]
                    work_queue.put((x_range, y_range))

            # create a queue to pass to workers to store the results
            result_queue = mp.Queue()
            procs = []

            # spawn workers
            for i in range(n_process):
                worker = Worker_median(work_queue, result_queue)
                procs.append(worker)
                worker.start()
                time.sleep(0.01)

            # collect the results off the queue
            for i in range(n_process):
                result_queue.get()

            for p in procs:
                p.join()

            print("Median subtraction took %.2f seconds." % (time.time() - tic))

        elif do_recipe=='convolve':

            kernel_sigma = 3.  # Kernel sigma in pixels
            kernel_size = np.full(2, round(4 * kernel_sigma / 2.) * 2 + 1, dtype=np.int)  # Kernel size in pixels

            print("Generating gaussian kernel")
            print('Kernel_size: ', kernel_size, ' - Kernel_sigma: ', kernel_sigma)
            kernel_data = np.asarray(
                Gaussian2DKernel(kernel_sigma, x_size=kernel_size[0], y_size=kernel_size[1],
                    mode='center'))  #mode='integrate'))
            kernel_size = kernel_data.shape

            shared_kernel_base = mp.Array(ctypes.c_float, kernel_data.size)
            shared_kernel = np.ctypeslib.as_array(shared_kernel_base.get_obj())
            shared_kernel = shared_kernel.reshape(kernel_size[0], kernel_size[1])
            shared_kernel[:] = kernel_data

            work_queue = mp.Queue()
            grid_n = np.asarray(grid_n)
            grid_mesh = np.ceil(im_size * 1. / grid_n).astype(int)
            grid_x = np.append(np.arange(0, im_size[0], grid_mesh[0]), im_size[0])
            grid_y = np.append(np.arange(0, im_size[1], grid_mesh[1]), im_size[1])

            for i in range(0, grid_x.size - 1):
                for j in range(0, grid_y.size - 1):
                    if work_queue.full():
                        print("Oh no! Queue is full after only %d iterations" % j)

                    x_range = [grid_x[i], grid_x[i + 1]]
                    y_range = [grid_y[j], grid_y[j + 1]]
                    work_queue.put((x_range, y_range))

            # create a queue to pass to workers to store the results
            result_queue = mp.Queue()
            procs = []

            # spawn workers
            for i in range(n_process):
                worker = Worker_convolve(work_queue, result_queue)
                procs.append(worker)
                worker.start()

            # collect the results off the queue
            for i in range(n_process):
                result_queue.get()

            for p in procs:
                p.join()

            print("Convolution took %.2f seconds." % (time.time() - tic))

        if os.path.isfile(im_file_out):
            print('Deleting existing image: ', im_file_out)
            os.remove(im_file_out)

        print('Writing image: ', im_file_out)
        pyfits.writeto(im_file_out, shared_nim, header=im_h)
