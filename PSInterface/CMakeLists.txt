#
# $Id$
#

# ROOT and pybindings didn't like each other
# two problems occur:
# 1) compilation error due to ambigious type_info types
#    in ROOT and boost:python
# 2) python crashes when loading the PSInterface wrapper
#    when ROOT is linke into libPSInterface
#    
# Here's a flag to remove the ROOT dependencies
# in order to make the pybindings work
# If you want pybindings set it to TRUE
set(DISABLE_ROOT FALSE)

i3_project(PSInterface)

if(DISABLE_ROOT)
  
  add_definitions(-DPSI_DISABLE_ROOT)

  i3_add_library(PSInterface 
    private/PSInterface/*.cxx private/PSInterface/*/*.cxx

    USE_TOOLS photonics ptd gsl boost log4cplus root

    USE_PROJECTS icetray phys-services dataclasses 
  )

  # only build pybindings if ROOT is disabled
  add_subdirectory(private/pybindings)

else(DISABLE_ROOT)

  i3_add_library(PSInterface 
    private/PSInterface/*.cxx private/PSInterface/*/*.cxx

    USE_TOOLS photonics ptd gsl boost log4cplus root

    USE_PROJECTS icetray phys-services dataclasses 

    ROOTCINT public/PSInterface/PSInterface.h 
    public/PSInterface/PSI_Coordinate.h 
    public/PSInterface/PSI_Dummy.h 
    public/PSInterface/PSI_FATInjector.h 
    public/PSInterface/PSI_PTD.h 
    public/PSInterface/PSI_Coordinate_PTD.h 
    public/PSInterface/PSI_Photonics.h 
    public/PSInterface/PSI_PhotonicsPhotorec.h 
    public/PSInterface/PSI_Coordinate_Photonics.h 
  )
endif(DISABLE_ROOT)


i3_test_executable(test 
  private/test/*.cxx
  USE_PROJECTS PSInterface
)
i3_test_scripts(resources/scripts/*.py)

