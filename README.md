ProcessNeuroImage
=================

Process modeling for functional and structural neuroimaging data

This project is based off the works by Andrew F. Hayes:
http://www.afhayes.com/

More specifically, he has written code to implement a large number of models using software packages like SPSS and SAS. The current package is implementing his ideas and models for use with neuroimaging data. Therefore, the current work is implementing his models, which are all outlined in the paper: http://www.afhayes.com/public/process2012.pdf.

Another key feature of the current work is that it is being written in MatLab for use on a computing cluster without the use of the distributed computing toolbox. The code uses the SGE and unix calls like "qsub." Without the use of a computing cluster these models may not be able to be estimated due to the intense computational time. The large computational time comes mainly from the fact that a large number of voxels (~100,000) are estimated in a typical brain and significance is tested using bootstrapping.