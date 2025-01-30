CC = gcc

SRCS = bloch_eul.c bloch_rk.c
OFILES = $(SRCS:.c=.o)

NAME = UM_blochsim.so

all: preclean $(NAME) oclean

$(NAME): $(OFILES)
	$(CC) -shared -o $(NAME) -fPIC $(SRCS)

oclean:
	rm -f $(OFILES)

preclean:
	rm -f $(NAME)
	rm -f $(OFILES)

clean: preclean