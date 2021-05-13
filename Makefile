# Compilers
CC = gcc

# Flags
EXTRA_CFLAGS=
CFLAGS = -std=c11 -Wall -Wextra $(EXTRA_CFLAGS)
LDFLAGS = -lm
CFLAGS_VECTOR = -O3 -march=native
# Binary file
TARGET = tiny_mc

# Files
C_SOURCES = tiny_mc.c wtime.c mtwister.c
C_OBJS = $(patsubst %.c, %.o, $(C_SOURCES))
C_SOURCES_VEC = wtime.c mtwister.c tiny_mc_vectorized.c
# Rules
all: $(TARGET)

$(TARGET): $(C_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

vector:
	$(CC) $(CFLAGS_VECTOR) -o tiny_mc_vectorized $(C_SOURCES_VEC) $(LDFLAGS)


clean:
	rm -f $(TARGET) *.o
