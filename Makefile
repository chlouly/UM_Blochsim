CC = gcc
CFLAGS = -Wall -Wextra
INCL = -Iincl

SRCS = src/bloch_eul.c src/bloch_rk4.c src/bloch_ljn.c
OFILES = $(SRCS:.c=.o)

NAME = UM_Blochsim.so

all: preclean $(NAME) oclean

%.o : %.c
	$(CC) $(INCL) -fPIC -c $< -o $@

$(NAME): $(OFILES)
	$(CC) -shared -o $(NAME) $(OFILES)

oclean:
	rm -f $(OFILES)

preclean:
	rm -f $(NAME)
	rm -f $(OFILES)

clean: preclean
