## This code was creating for showing the use of shared memory and multiprocess in Python3

*Implemented recipes:*
-median
-convolve


*Flags*
-ncpu: Number of CPUs
-ncode: Number of cores
-nx: Number of bins along X for processing image
-ny: Number of bins along Y for processing image
-fft: Apply the FFT convolution
-fft_interpol: Interpolate NaN values when doing FFT convolution

*Examples*

python3 process_image.py --recipe convolve -f --ncore 1 --nx 1 --ny 1

python3 process_image.py --recipe convolve -f --ncore 2 --nx 2 --ny 2

python3 process_image.py --recipe convolve -f --ncore 2 --nx 2 --ny 2 --fft --fft_interpol

python3 process_image.py --recipe convolve -f --ncore 2 --nx 2 --ny 2 --fft

Author: Roberto Pablo Muñoz
Email: rmunoz@uc.cl