# Compilers
CC = gcc

# Flags
EXTRA_CFLAGS=
CFLAGS = -std=c11 -Wall -Wextra $(EXTRA_CFLAGS)
LDFLAGS = -lm
CFLAGS_VECTOR = -O3 -march=native $(EXTRA_CFLAGS)
CFLAGS_OMP_LINEAL = -O3 -march=native -flto -fopenmp $(EXTRA_CFLAGS)
CFLAGS_OMP_VECTOR = -O3 -ffast-math -flto -march=native -fopenmp -fopenmp-simd $(EXTRA_CFLAGS)

# Binary file
TARGET = tiny_mc

# Files
C_SOURCES = tiny_mc.c wtime.c mtwister.c
C_OBJS = $(patsubst %.c, %.o, $(C_SOURCES))
C_SOURCES_VEC = wtime.c mtwister.c tiny_mc_vectorized.c
C_SOURCES_OMP = wtime.c mtwister.c tiny_mc_omp.c
C_SOURCES_VEC_OMP = wtime.c mtwister.c tiny_mc_vectorized_omp.c
# Rules
all: $(TARGET)

$(TARGET): $(C_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

vector:
	$(CC) $(CFLAGS_VECTOR) -o tiny_mc_vectorized $(C_SOURCES_VEC) $(LDFLAGS)

omp_vector:
	$(CC) $(CFLAGS_OMP_VECTOR) -o tiny_mc_vectorized_omp $(C_SOURCES_VEC_OMP) $(LDFLAGS)

omp_lineal:
	$(CC) $(CFLAGS_OMP_LINEAL) -o tiny_mc_omp $(C_SOURCES_OMP) $(LDFLAGS)

clean:
	rm -f $(TARGET) *.o
