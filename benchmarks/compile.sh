./configure --enable-single --enable-threads \
	--enable-sse --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-avx-128-fma \
	--enable-generic-simd128  --enable-generic-simd256

# under benchfft-3.1
# $: ./configure ...
# $: make

# under benchfft-3.1/benchees/fftw3
# $: make benchmark

# Raw Benchmark Data
# name-of-code transform-type transform-size mflops time setup-time

# Benchmark type can be configured in scripts/benchmark and benchees/fftw3/Makefile