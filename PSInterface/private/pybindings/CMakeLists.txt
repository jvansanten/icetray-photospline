
if(BUILD_PYBINDINGS)

  i3_add_pybindings(PSInterface
    PSInterface.cxx
    USE_TOOLS boost python log4cplus
    USE_PROJECTS PSInterface icetray
    )
endif(BUILD_PYBINDINGS)

