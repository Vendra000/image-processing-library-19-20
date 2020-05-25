test_bmp: bmp.o test_bmp.c
	gcc test_bmp.c bmp.o -o test_bmp -Wall -lm

bmp.o: bmp.c
	gcc bmp.c -o bmp.o -c -Wall -lm

ip_lib.o: bmp.c ip_lib.c
	gcc -c bmp.c ip_lib.c -Wall -lm

test: bmp.o ip_lib.o main_iplib.c
	gcc bmp.o ip_lib.o main_iplib.c -Wall --ansi --pedantic -lm -g3 -o test -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra
