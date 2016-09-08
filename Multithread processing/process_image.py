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
            kernel_shape = np.asarray(shared_kernel.shape)

            dx=(x_range[1]-x_range[0])
            dy=(y_range[1]-y_range[0])
            im_padx = np.full(2, np.ceil( (kernel_shape[1] - 1) / 2), dtype=int)
            im_pady = np.full(2, np.ceil( (kernel_shape[0] - 1) / 2), dtype=int)

            if x_range[0]==0:
                im_padx[0]=0
            if x_range[1]==shared_im_shape[1]:
                im_padx[1]=0

            if y_range[0]==0:
                im_pady[0]=0
            if y_range[1]==shared_im_shape[0]:
                im_pady[1]=0

            if do_debug:
                print('IMAGE SHAPE: ', shared_im_shape)
                print('IMAGE RANGE - x_range: ', x_range, ' - y_range: ', y_range)
                print('dx: ', dx, ' - im_padx: ', im_padx)
                print('dy: ', dy, ' - im_pady: ', im_pady)
                print('CONVOLVE_FFT - x_range: ', [x_range[0] - im_padx[0], x_range[1]+im_padx[1]], ' - y_range: ', [y_range[0] - im_pady[0], y_range[1] + im_pady[1]])
                print('CROP -  x_range: ', [im_padx[0], dx+im_padx[0]], ' - y_range: ', [im_pady[0], dy+im_pady[0]])
                print('')

            shared_nim[y_range[0]:y_range[1], x_range[0]:x_range[1]] = convolve(
                shared_im[y_range[0] - im_pady[0]:y_range[1] + im_pady[1], x_range[0] - im_padx[0]:x_range[1] + im_padx[1]],
                shared_kernel, normalize_kernel=True)[im_pady[0]:dy + im_pady[0], im_padx[0]:dx + im_padx[0]]

            #print('Worker_convolve done: ', mp.current_process())
            self.result_queue.put(id)

class Worker_convolve_fft(mp.Process):
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

            kernel_shape = np.asarray(shared_kernel.shape)

            dx=(x_range[1]-x_range[0])
            dy=(y_range[1]-y_range[0])
            im_padx = np.full(2, np.ceil( (kernel_shape[1] - 1) / 2), dtype=int)
            im_pady = np.full(2, np.ceil( (kernel_shape[0] - 1) / 2), dtype=int)

            if do_padding:
                if x_range[0]==0 and x_range[1]==shared_im_shape[1]:
                    im_padx[0]=0
                    im_padx[1]=0
                else: 
                    if x_range[0]==0:
                        im_padx[0]=0
                    if x_range[1]==shared_im_shape[0]:
                        im_padx[1]=0

                    if x_range[0]>0:
                        im_padx[0]= int(2**(np.ceil(np.log2(dx+im_padx[0]+im_padx[1]))) - (dx+im_padx[1]))
                    if x_range[1]<shared_im_shape[0]:
                        im_padx[1]= int(2**(np.ceil(np.log2(dx+im_padx[0]+im_padx[1]))) - (dx+im_padx[0]))

                if y_range[0]==0 and y_range[1]==shared_im_shape[0]:
                    im_pady[0]=0
                    im_pady[1]=0
                else:
                    if y_range[0]==0:
                        im_pady[0]=0
                    if y_range[1]==shared_im_shape[0]:
                        im_pady[1]=0

                    if y_range[0]>0:
                        im_pady[0]= int(2**(np.ceil(np.log2(dy+im_pady[0]+im_pady[1]))) - (dy+im_pady[1])) 
                    if y_range[1]<shared_im_shape[0]:
                        im_pady[1]= int(2**(np.ceil(np.log2(dy+im_pady[0]+im_pady[1]))) - (dy+im_pady[0]))   

            else:
                if x_range[0]==0:
                    im_padx[0]=0
                if x_range[1]==shared_im_shape[1]:
                    im_padx[1]=0

                if y_range[0]==0:
                    im_pady[0]=0
                if y_range[1]==shared_im_shape[0]:
                    im_pady[1]=0

            if do_debug:
                print('IMAGE SHAPE: ', shared_im_shape)
                print('IMAGE RANGE - x_range: ', x_range, ' - y_range: ', y_range)
                print('dx: ', dx, ' - im_padx: ', im_padx)
                print('dy: ', dy, ' - im_pady: ', im_pady)
                print('CONVOLVE_FFT - x_range: ', [x_range[0] - im_padx[0], x_range[1]+im_padx[1]], ' - y_range: ', [y_range[0] - im_pady[0], y_range[1] + im_pady[1]])
                print('CROP -  x_range: ', [im_padx[0], dx+im_padx[0]], ' - y_range: ', [im_pady[0], dy+im_pady[0]])
                print('FFT_INTERPOL: ', do_fft_interpol)
                print('')

            shared_nim[y_range[0]:y_range[1], x_range[0]:x_range[1]] = convolve_fft(
                shared_im[y_range[0] - im_pady[0]:y_range[1] + im_pady[1], x_range[0] - im_padx[0]:x_range[1] + im_padx[1]],
                shared_kernel, normalize_kernel=True, interpolate_nan=do_fft_interpol)[im_pady[0]:dy + im_pady[0],
                im_padx[0]:dx + im_padx[0]]

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

            shared_nim[y_range[0]:y_range[1], x_range[0]:x_range[1]] = \
                shared_im[y_range[0]:y_range[1], x_range[0]:x_range[1]] - \
                np.nanmedian(shared_im[y_range[0]:y_range[1], x_range[0]:x_range[1]])

            #print('Worker_median done - id: ', mp.current_process().name)
            self.result_queue.put(id)


if __name__ == "__main__":

    argv=sys.argv[1:]
    im_file_in='data/VCC1010.MASKED.fits'
    do_recipe='median'
    do_overwrite=False
    do_fft=False
    do_fft_interpol=False
    do_padding=False
    do_debug=False
    n_cpu = 1
    n_core = 4
    n_x=1
    n_y=6
    hyper_thread=1 # If you have hyper threading, then set the value to 2, otherwise set to 1

    kernel_sigma = 3.  # Kernel sigma in pixels
    kernel_shape = np.full(2, round(4 * kernel_sigma / 2.) * 2 + 1, dtype=np.int)  # Kernel size in pixels

    try:
        opts, args = getopt.getopt(argv, "hfri:o:", ["recipe=", "ifile=", "ofile=","ncpu=","ncore=","nx=","ny=","fft","padding","debug","fft_interpol"])
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
        elif opt == '--fft':
            do_fft=True
        elif opt == '--fft_interpol':
            do_fft_interpol=True
        elif opt == '--padding':
            do_padding=True
        elif opt == '--debug':
            do_debug=True

    n_process = np.min( [n_cpu * n_core * hyper_thread, n_x*n_y] )
    grid_n = [n_y, n_x]  # Number of bins to divide image along [Y,X]
    im_file_out = 'data/VCC1010.MASKED.'+do_recipe.upper()+'.fits'

    print('Reading image: ', im_file_in)
    hdulist = pyfits.open(im_file_in)
    im_data = hdulist[0].data
    im_h = hdulist[0].header
    hdulist.close()
    im_shape = np.asarray(im_data.shape)
    im_size = int(np.prod(im_shape))
    print('Image size: ', im_size)

    if do_padding:
        print('Adding padding image to boost FFT Convolution')
        shared_im_shape=2**np.ceil(np.log2(im_shape))
        if n_y>1:
            shared_im_shape[0]=im_shape[0]
        if n_x>1:
            shared_im_shape[1]=im_shape[1]
    else:
        shared_im_shape=im_shape

    shared_im_size=int(np.prod(shared_im_shape))

    shared_im_base = mp.Array(ctypes.c_float, shared_im_size)
    shared_im = np.ctypeslib.as_array(shared_im_base.get_obj())
    shared_im = shared_im.reshape(shared_im_shape[0], shared_im_shape[1])
    shared_im[:] = np.nan
    shared_im[0:im_shape[0],0:im_shape[1]] = im_data

    shared_nim_base = mp.Array(ctypes.c_float, shared_im_size)
    shared_nim = np.ctypeslib.as_array(shared_nim_base.get_obj())
    shared_nim = shared_nim.reshape(shared_im_shape[0], shared_im_shape[1])
    shared_nim[:] = np.nan

    if os.path.isfile(im_file_out) is False or do_overwrite is True:

        tic = time.time()
        if do_recipe=='median':

            work_queue = mp.Queue()
            grid_n = np.asarray(grid_n)
            grid_mesh = np.ceil(shared_im_shape * 1. / grid_n).astype(int)
            grid_x = np.append(np.arange(0, shared_im_shape[1], grid_mesh[1]), shared_im_shape[1])
            grid_y = np.append(np.arange(0, shared_im_shape[0], grid_mesh[0]), shared_im_shape[0])

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

            if do_padding:
                print('Removing padding')
                print('shared_nim BEFORE: ', shared_nim.shape)
                shared_nim=shared_nim[0:im_shape[0],0:im_shape[1]]
                print('shared_nim AFTER: ', shared_nim.shape)

            print("Median subtraction took %.2f seconds." % (time.time() - tic))

        elif do_recipe=='convolve':          

            print("Generating gaussian kernel")
            print('Kernel_size: ', kernel_shape, ' - Kernel_sigma: ', kernel_sigma)
            kernel_data = np.asarray(
                Gaussian2DKernel(kernel_sigma, x_size=kernel_shape[0], y_size=kernel_shape[1],
                    mode='center'))  #mode='integrate'))
            kernel_shape = kernel_data.shape

            shared_kernel_base = mp.Array(ctypes.c_float, kernel_data.size)
            shared_kernel = np.ctypeslib.as_array(shared_kernel_base.get_obj())
            shared_kernel = shared_kernel.reshape(kernel_shape[0], kernel_shape[1])
            shared_kernel[:] = kernel_data

            work_queue = mp.Queue()
            grid_n = np.asarray(grid_n)
            grid_mesh = np.ceil(shared_im_shape * 1. / grid_n).astype(int)
            grid_x = np.append(np.arange(0, shared_im_shape[1], grid_mesh[1]), shared_im_shape[1])
            grid_y = np.append(np.arange(0, shared_im_shape[0], grid_mesh[0]), shared_im_shape[0])

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
                if do_fft:
                    worker = Worker_convolve_fft(work_queue, result_queue)
                else:
                    worker = Worker_convolve(work_queue, result_queue)
                procs.append(worker)
                worker.start()

            # collect the results off the queue
            for i in range(n_process):
                result_queue.get()

            for p in procs:
                p.join()

            if do_padding:
                print('Removing padding')
                print('shared_nim BEFORE: ', shared_nim.shape)
                shared_nim=shared_nim[0:im_shape[0],0:im_shape[1]]
                print('shared_nim AFTER: ', shared_nim.shape)
            
            print("Convolution took %.2f seconds." % (time.time() - tic))

        if os.path.isfile(im_file_out):
            print('Deleting existing image: ', im_file_out)
            os.remove(im_file_out)

        print('Writing image: ', im_file_out)
        pyfits.writeto(im_file_out, shared_nim, header=im_h)
