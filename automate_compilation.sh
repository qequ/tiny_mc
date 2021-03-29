cat parameters.txt | while read p; do for ((i=1; i<=30; i=i+1)); do
make clean && 
make EXTRA_CFLAGS="${p}" CPP_FLAGS="-DPHOTONS=$((1 + $RANDOM % 32768))" && 
./tiny_mc; done; done