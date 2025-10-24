CC = gcc
CFLAGS = -O2 -std=c11 -Wall -pthread -Iinclude

SRC = src/main.c src/gauss.c
OBJ = $(SRC:.c=.o)
TARGET = gauss_mt

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
