gcc -c compute.c
gcc -shared -o compute.dll compute.o
pyinstaller --onefile --noconsole --add-data compute.dll;. main.py