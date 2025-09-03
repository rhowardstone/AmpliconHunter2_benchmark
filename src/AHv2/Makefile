CC = gcc
CFLAGS = -Wall -O3 -march=native -pthread -ffast-math
LDFLAGS = -pthread -lm
TARGET = amplicon_hunter

# Source files
SRCS = amplicon_hunter.c
OBJS = $(SRCS:.c=.o)

# Default target
all: $(TARGET)

# Build the executable
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Debug build
debug: CFLAGS = -Wall -g -O0 -pthread -DDEBUG -fsanitize=thread
debug: LDFLAGS = -pthread -lm -fsanitize=thread
debug: clean $(TARGET)

# Profile-guided optimization
pgo: clean
	$(CC) $(CFLAGS) -fprofile-generate -o $(TARGET) $(SRCS) $(LDFLAGS)
	@echo "Run typical workloads, then 'make pgo-use'"

pgo-use:
	$(CC) $(CFLAGS) -fprofile-use -o $(TARGET) $(SRCS) $(LDFLAGS)

# Ultra-optimized build
ultra: CFLAGS = -Wall -Ofast -march=native -pthread -ffast-math -funroll-loops -flto
ultra: LDFLAGS = -pthread -lm -flto
ultra: clean $(TARGET)

# Clean
clean:
	rm -f $(TARGET) $(OBJS) *.gcda *.gcno

# Install
install: $(TARGET)
	install -m 755 $(TARGET) /usr/local/bin/

# Test compress
test-compress:
	./$(TARGET) compress --input-dir test/G100 --output test/G100_compressed --threads 8

# Test run
test-run:
	./$(TARGET) run --input test/file_list.txt --primers test/primers.txt \
		--output test_output.fa --threads 8 --mismatches 2 --clamp 3 \
		--min-length 1000 --max-length 2000

.PHONY: all clean debug pgo pgo-use ultra install test-compress test-run