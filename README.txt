# Compiler optimization using AI

Project by Mairaj Khalid, student ID: 20019274
For King's College London Computer science with intelligent systems BS.c final year project 2022-2023

------------------------------

## Project files/folders:

------------------------------

## Items included are:

- main.cpp // Main source file.
- include // Folder with any needed include files for compilation.
- compiler_optimization_flags_information.csv // CSV file with required items for program.exe during execution.
- compiler_optimization_flags_information_not_used.txt // Text file with extra information.
- compilation_commands_for_test_cases.txt // Text file with required items for program.exe during execution.
- test_programs // Source files for testing compilation that are within two layers of folders.
- README.txt // README of the project and current file.

------------------------------

## Programs the project was built using and how the project was built:

------------------------------

- 1: Processor used was Intel 8th gen corei7 4 core processor with Intel UHD graphics.
- 2: The OS environment built upon was Windows 10.
- 3: Built using IDE Code::blocks version 20.03.
- 4: The compiler used was MinGW-W64 gcc version 8.1, provided by Code::blocks. Link (should only work for Windows 10 OS):

https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/8.1.0/

- 5: g++ C++11 was used for compilation.
- 6: The compiler flags used for optimization are based on the gcc version 8.1, and are inside compiler_optimization_flags_information.csv, while those not
used are inside compiler_optimization_flags_information_not_used.txt (See "Extra file information" section on how to change them).
Link (For gcc flags summary and optimization flags full explanation):

https://gcc.gnu.org/onlinedocs/gcc-8.1.0/gcc/Option-Summary.html#Option-Summary
https://gcc.gnu.org/onlinedocs/gcc-8.1.0/gcc/Optimize-Options.html#Optimize-Options

- 7: Test cases used are from MiBench 1.0. Link:

https://vhosts.eecs.umich.edu/mibench/source.html

- 8: The current program and source code are tuned for the test programs. This can be easily changed by adjusting constant values inside the main.cpp
file for use-case (See "Changing variable values to modify algorithm").
- 9: There is ability to plot graphs, but this is disabled and not required in current executable and source code (See "Enabling plotting" section
to enable it).
- 10: The ability to plot graphs is by using matplotlib and numpy from Python, but for C++ matplotlibcpp.h file (inside include folder) was used. This file
acted as a bridge between Python and C++ to use matplotlib and numpy.
- 11: Python version 3.8.10 was used for its matplotlib and numpy for plotting graphs only.

## Note: Python and matplotlibcpp are NOT needed for current state of project, unless enabling plotting.

------------------------------

## Installation requirements:

------------------------------

- 1: Windows 10 OS.
- 2: gcc compiler version 8.1 or higher with g++ with C++11 available.

------------------------------

## Compiling source code:

------------------------------

- 1: If completed "Installation requirements" section, execute command while within source folder:

g++ -O3 -std=c++11 main.cpp -o program

------------------------------

## Running program:

------------------------------

## Inputs:

- 1: If completed "Installation requirements" and "Compiling source code" sections, execute program.exe.
- 2: Input compilation command you wish to optimize for when prompted (See "Running with test cases section" for getting pre-set commands).
- 3: Input location of executable created on compilation when prompted (See "Running with test cases section" for getting pre-set locations).

## Note: Executing program.exe should work on Windows 10 (if not then see "Compiling source code" section).

Outputs:

- 1: Program will entirely output on command line (See "Enabling plotting" section if graphs are required).
- 2: Outputs are following, and separated by dashed lines:
- 2.1: Inputs prompts.
- 2.2: Notifications such as initialization start with average run-time, main loop start, and final outputs reached.
- 2.3: Current pass that ended on main loop.
- 2.4: Standard deviation, number of re-runs of Chromosomes, and average run-time after execution on current loop.
- 2.5: Final outputs such as base run-time, flag to compare against, final flags chosen, etc.
- 3: Set errors if part 2 and 3 of "Inputs" section are incorrect.

------------------------------

## Running with test cases:

------------------------------

- 1: Choose whichever test case you wish to go with inside compilation_commands_for_test_cases.txt.
- 2: Note the that a compilation command has a matching output executable file path noted lower.
- 3: Based on your selection use those as inputs for the "Inputs" section of "Running program.exe" section.

#Note: There are other test cases available to test against inside test_programs folder, but their compilation command and matching output file paths are
not set up (though current ones should be enough).

------------------------------

## Changing variable values to modify algorithm:

------------------------------

- 1: These changes to be made in main.cpp.
- 2: All changes to affect the algorithm are the set up constant values with bold letters. they have comments on them to describe what values are allowed.
- 3: Recompiling will be necessary to see these changes in the executable (See "Compiling source code" section).

------------------------------

## Extra file information:

------------------------------

- 1: compiler_optimization_flags_information.csv: This contains a list of flags and the values they can take. Flags that can be turned on and off only have
just themselves in a single line, while flags that can take multiple values have themselves ending with equal sign "=" and then the values they can take
separated by commas "," on a single line. Just simply add or remove flags based on how they are setup.
- 2: compiler_optimization_flags_information_not_used.txt: This contains a list of flags that the compiler will not be using with explanations. This file
is only for the user, and will not be used by the program.

------------------------------

## Enabling plotting:

------------------------------

- 1: Install python version 3.X (preferred 3.8.10) on your system, and should be available in PATH environment.
- 2: Uncomment inside main.cpp file all code under the heading of "Plotting." comment, and save.
- 3: Now recompile program.exe wit command:

g++ -O3 -I<Path to Python include folder> -I<Path to numpy include folder> -I./include
-L<Path to Python libs folder> -lpython<Python version number (example: 38)> main.cpp -o program

- 4: Now simply execute program.exe as in "Running program" section.
