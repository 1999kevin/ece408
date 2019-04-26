# ece408
(This is under Windows 10 system)
Please follow the steps below:

1) Generate the input data
python export.py

2) Clean the compiled files
make clean

3) Compile the files
make

4) Do the computation
compute

*Note: The test_import.c is just used for debugging the input, if you want to test it, please enter the command below:
gcc -o test_import test_import.c import.c
test_import

##Important##
Please do not delete the "data" folder. It contains the input and output .csv files.
You can delete those files manually.
