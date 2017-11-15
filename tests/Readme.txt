================================================================================
=========================== potfit testing framework ===========================
================================================================================

This directory contains simple scripts to test the potfit runtime
behaviour for many different configurations.

--------------------------------------------------------------------------------

A. SYSTEM REQUIREMENTS

To run the tests the following software needs to be installed:

1. Python3 runtime with the pytest module
2. clang compiler version 3.8 or later

--------------------------------------------------------------------------------

B. INTRODUCTION

The tests are split up based on the potential model and interaction used in
the respective directory. Analytic pair potentials are in apot/pair.
Additional features like stress or evo are then in another subdirectory,
e.g. apot/pair/stress.

All tests which are related to the same feature go into the same file
which needs to start with 'test_' and continue with the model and interaction,
like 'test_apot_pair.py'.

A test can be defined as a regular python function starting with 'test_'.
The argument 'potfit' will be a python class provided by the pytest framework
and represent the potfit binary corresponding to the current directory. Check
the conftest.py file for details.

--------------------------------------------------------------------------------

C. EXAMPLE

  def test_apot_pair_sc(potfit):
      # Create a standard parameter file
      potfit.create_param_file()

      # Call makeapot to generate a basic potential file
      potfit.call_makeapot('startpot', '-n 1 -i pair -f eopp_sc')

      # Create a simple config file
      potfit.create_config_file(repeat_cell=3, seed=42)

      # Run the test
      potfit.run()

      # Examine the results
      assert potfit.has_no_error()

--------------------------------------------------------------------------------

D. API DOCUMENTATION

The potfit class has the following functions for setting up tests:

## create_param_file(keyword=value, ...)

  Creates a potfit parameter file name 'param_file' in the testing directory
  with the following default values:

    ntypes 1
    startpot startpot
    endpot endpot
    config config
    tempfile tempfile
    output_prefix output

  All values can be adjusted by passing them as named arguments:

    create_param_file(config='conf_file', opt=1, eng_weight=3.14159)

## create_potential_file(string)

  Creates a file called 'startpot' with the content of the string argument.

## create_config_file()

  Creates a rather simple config file with atoms randomly being displaced
  from their positions. It accepts the following named arguments:

    repeat_cell
        number of unit cells per direction (default 2)
    seed
        random number generator seed for displacements
    ntypes
        number of different atom types to generate (default 1)
    useforce
        set the useforce flag
    energy
        cohesive energy for the configuration (default -1.0)
    stress [bool]
        write some random stress value (default Off)

## run(param_file='param_file')

  Run potfit binary with the provided parameter file.
  Once this function returns the following class members will be available:

    potfit.stdout
        contains the potfit output to stdout
    potfit.stderr
        contains the potfit output to stderr
    potfit.startpot
        contains the content of the startpot file
    potfit.endpot
        contains the content of the endpot file
    potfit.force
        contains the content of the force output file
    potfit.energy
        contains the content of the energy output file
    potfit.stress
        contains the content of the stress output file
    potfit.punish
        contains the content of the punishments file
    potfit.rho_loc
        contains the content of the local electron density file

  If any of these files does not exist the content of the variables
  will be an empty string.

## has_no_error()

  Convenience function for checking if no error occurred:

    assert potfit.has_no_error()

  Checks for a return value of 0 and an empty stderr output.

## has_no_warning()

  Convenience function for checking if an warning was emitted:

    assert potfit.has_no_warning()

  Checks if '[WARNING]' is present in stdout.

## has_error()

  Convenience function for checking if an error occurred:

    assert potfit.has_error()

  Checks if the return value is != 0 and '[ERROR]' is present in stderr.

## create_file(filename)

  If the provided functions for creating input files are not enough this
  function can be used to create custom input files. Please note that
  this function will return an open file which needs to be closed before
  calling potfit.run().

    f = potfit.create_file('my_pot')
    f.write('....\n')
    f.close()
