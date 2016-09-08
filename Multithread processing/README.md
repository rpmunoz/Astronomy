** How to run the code

python3 process_image.py --recipe convolve -f --ncore 1 --nx 1 --ny 1

python3 process_image.py --recipe convolve -f --ncore 2 --nx 2 --ny 2

python3 process_image.py --recipe convolve -f --ncore 2 --nx 2 --ny 2 --fft --fft_interpol

python3 process_image.py --recipe convolve -f --ncore 2 --nx 2 --ny 2 --fft