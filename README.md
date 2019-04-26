For Windows 10 system:
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


For Mac, we can use Xcode to compile and run it:

1)open Xcode and create a new Xcode project(macOS template: Command Line Tool)

2)add all files into your working directory except test_import.c and delete the original main.c

3)click Run.


##Important##
Please do not delete the "data" folder. It contains the input and output .csv files.
You can delete those files manually.
