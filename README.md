# Welcome to maserVlbi
a set of classes and tools to handle the VLBI maser data.

Set consists of:

### CLASSES ###

- maserVlbi
    -  spotClass
    -  spectrumClass
    -  table of cloudletClass
        - spotClass
        - calculated properties


### TOOLS ###
- asciiToHDF5Converter
- cloudletFinder

### USAGE ###
asciiToHDF5Converter -conf exampleConfig.ini -o file.hdf5
cloudletFinder file.hdf5
