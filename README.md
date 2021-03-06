# φasma (gr-phasma)

### Introduction
Φasma is a GNURadio module that aims at the contribution to the research in the field of Signal Detection and Automatic Modulation Classification by providing various signal processing and machine learning blocks.

### Installation
To install gr-phasma, the following dependencies are required:
- [GNU Radio](https://github.com/gnuradio/gnuradio)
- [Armadillo](http://arma.sourceforge.net/download.html)
- [OpenCV](http://opencv.org/)
- [JsonCpp](https://github.com/open-source-parsers/jsoncpp)


Install by the following shell commands:
```
git clone git@github.com:ctriant/gr-phasma.git
cd gr-phasma
mkdir build
cd build
cmake ..
make -j4
sudo make install
```

### Usage
For detailed information about each block, visit the Doxygen HTML pages.

### Brief block description
Below there is a brief description of the basic blocks provided by gr-phasma

##### Eigenvalue Signal Detector

##### Signal Extractor

##### OpenCV Random Forest Modeler

##### OpenCV Classifier


### License
This software is Copyright © 2016 Kostis Triantafyllakis. It is free software, and is released under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) License.

### Contact
Maintainer of this module:

Kostis Triantafyllakis<br/>
ctriant[at]csd.uoc.gr
