#
# This script works only with Python 3
#
import os
import re

import create_initial_conditions
import Test

MaxParticlesInSequentialUpscalingStudies = 200


def step1(compiler):
  test = Test.Test( compiler, "step-1" )
  test.cleanall()

  arguments = "CXXFLAGS=\"--std=c++0x\""
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  _,arguments = create_initial_conditions.create_grid_setup( 2,1,1, 1, 1, "no-noise" , 0.1, 0.1, 1 )
  print( "Run code with " + arguments )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()

  search_pattern = "([-+]?\d*\.\d+|\d+), *([-+]?\d*\.\d+|\d+), *([-+]?\d*\.\d+|\d+)"
  result = test.search_for_pattern_in_output(search_pattern)
  if result!="":
    print( "Got " + result + " as last output which I interprete to be a particle position ... ok (though that does not mean that the data is correct; that's something I don't validate here)" )
  else:
    print( "Last line of output should still be the one I used in the template ... failed" )
    exit()

def step2(compiler):
  test = Test.Test( compiler, "step-2" )

  arguments = "CXXFLAGS=\"-O3 -fopenmp --std=c++0x\""
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  _,arguments = create_initial_conditions.create_grid_setup( 2,1,1, 1, 1, "no_noise" , 0.1, 0.1, 1 )
  print( "Run code with " + arguments )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok (but this does not mean that the outcome is correct)" )
  else:
    print( "Run failed: " + test.last_output )
    exit()

def step3(compiler):
  test = Test.Test( compiler, "step-3" )
  test.cleanall()

  arguments = "CXXFLAGS=\"-fno-tree-vectorize\""
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  particles_counts = [11,11,11]
  _,arguments = create_initial_conditions.create_grid_setup( particles_counts[0],particles_counts[1],particles_counts[2], 1, 1, "random" , 0.1, 0.0001, 10 )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()
  no_vec_time = test.runtime
  test.cleanall()

  arguments = "CXXFLAGS=\"-O3 -fopenmp --std=c++0x\""
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  particles_counts = [11,11,11]
  _,arguments = create_initial_conditions.create_grid_setup( particles_counts[0],particles_counts[1],particles_counts[2], 1, 1, "random" , 0.1, 0.0001, 10 )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()
  with_vec_time = test.runtime

  if no_vec_time<=with_vec_time:
    print( "Code is slower with vectorisation, so you might want to tune it ... check" )
    exit()
  else:
    print( "Code is already faster by a factor of " + str(no_vec_time/with_vec_time) + " through vectorisation but you might want to tune it further ... ok" )


def step4(compiler):
  test = Test.Test( compiler, "step-4" )

  arguments = "CXXFLAGS=\"-O3 --std=c++0x -fopenmp\""
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  particles_counts = [11,11,11]
  _,arguments = create_initial_conditions.create_grid_setup( particles_counts[0],particles_counts[1],particles_counts[2], 1, 1, "random" , 0.1, 0.0001, 10 )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )

  print( "Test for one core" )
  result = test.run( arguments, {"OMP_NUM_THREADS": "1"} )
  if result==1:
    print( "Run terminated successfully after " + str(test.runtime) + "s ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()
  serial_runtime = test.runtime

  for p in range(2,28,2):
    print( "Test for " + str(p) + " cores" )
    result = test.run( arguments, {"OMP_NUM_THREADS": str(p)} )
    speedup = serial_runtime/test.runtime
    if result==1 and speedup>p*0.8:
      print( "Run terminated successfully after " + str(test.runtime) + "s, i.e. with speedup of " + str(speedup) + " ... ok  (but you might want to tune it further)" )
    elif result==1:
      print( "Run terminated successfully after " + str(test.runtime) + "s, i.e. without expected speedup ... check runtimes" )
    else:
      print( "Run failed: " + test.last_output )
      exit()




if __name__=="__main__":
	print("""
 This is a small Python script to check that the format of a submission is
 valid. It does not check for correctness, although it gives some clues on
 performance. If you use the script on Hamilton, which is the machine
 we will use to mark the submission, you need to load appropriate modules:

 module load python/3.9.9 intel/2021.4;

 We recommended that you use a compute node rather than the login nodes. You
 can do this with the command:

 salloc -N 1 -p test.q python3 validate.py

 Disclaimer: This is only a sanity check, and we'll run further tests on
 correctness and scalability to mark the coursework. But the script ensures
 that your submission is formally correct and performs some basic checks.

""")
	compiler="gcc"
	try:
		print("Checking step 1")
		step1(compiler)
		print("Checking step 2")
		step2(compiler)
		print("Checking step 3")
		step3(compiler)
		print("Checking step 4")
		step4(compiler)

	except BaseException as err:
   		print(f"Unexpected {err=}, {type(err)=}")
   		raise
