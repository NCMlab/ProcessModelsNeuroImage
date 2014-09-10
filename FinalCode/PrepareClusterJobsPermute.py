import os
import sys
import glob
import datetime
import numpy as np

def main():
    if sys.argv is None:
        print 'Usage: PrepareClusterJobsPermute.py <BaseDir> <NPermute> <NJobSplit>'
        sys.exit(0)
    
    BaseDir = sys.argv[1]
    NPermute = sys.argv[2]
    NJobSplit =  sys.argv[3]
    
    if os.path.exists(BaseDir):
        DataDir=os.path.join(BaseDir,'data')
        if os.path.exists(DataDir):
            JobDir = os.path.join(BaseDir,'JobFiles')
            JobOutDir = os.path.join(BaseDir,'JobOutput')
            if not os.path.exists(JobDir):
                os.mkdir(JobDir)
            # Data file
            DataFile = glob.glob(os.path.join(DataDir,'ModelInfo.mat'))

            # create the file that will be a list of qsub calls
            QSubFile=os.path.join(JobDir,('QSubFile.sh'))
            Qfid=open(QSubFile,'write')
            # Number of permutes per job
            NPermPerJob = np.ceil(float(NPermute)/NJobSplit)
            for i in range(0,NJobSplit,1):

                JobFileName = os.path.join(JobDir,('job_%s.sh'%(i)))
                ParamFileName = os.path.join(JobDir,('param_%s.in'%(i)))
                # create the parameter files
                # remove the .mat from the end of the filename
                
                fid=open(ParamFileName,'write')
                fid.write('InputData=%s'%(DataFile[:-4]))
                fid.write('Count=%s'%(i))
                fid.write('NPermPerJob=%s'%(NPermPerJob))
                
                fid.close()
                # Create the job files
                fid=open(JobFileName,'write')
                fid.write('#!/bin/bash\n')
                fid.write('#PBS -S /bin/bash\n\n')
                fid.write('#Choose MCR directory\n')
                fid.write('MCR=/usr/local/matlab-mcr/v80\n\n')
                fid.write('echo "Running on host: `hostname`"\n')
                fid.write('cd $PBS_O_WORKDIR\n')
                fid.write('echo "Current working directory is `pwd`"\n\n')
                fid.write('echo "Starting run at: `date`"\n\n')
                fid.write('/home/steffejr/scripts/ProcessModelsNeuroImage/CompiledCode/run_RunCompiledPermutation.sh $MCR %s\n'%ParamFileName)
                fid.write('echo "Job finished at: `date`"\n\n')
                fid.close()
                # For this to work the input to this compiled program has to be a parameter file
                Qfid.write('qsub -e %s -o %s -l walltime=72:00:00:00,mem=3gb %s\n'%(JobOutDir,JobOutDir,JobFileName))
                
                #Str='qsub -l walltime=72:00:00:00,mem=6gb %s'%JobFileName 
                #os.system(Str)
            Qfid.close()
if __name__=='__main__':
    main()        

  
        
        
