import os
import sys
import glob
import datetime
#class FSExtract():
# def Extract(StudyName):
def main():
    if sys.argv is None:
        print 'Usage: PrepareClusterJobs.py <BaseDir>'
        sys.exit(0)
    
    BaseDir = sys.argv[1]
    if os.path.exists(BaseDir):
        DataDir=os.path.join(BaseDir,'data')
        if os.path.exists(DataDir):
            JobDir=os.path.join(BaseDir,'JobFiles')
            if not os.path.exists(JobDir):
                os.mkdir(JobDir)
            # How many data files
            DataFiles=glob.glob(os.path.join(DataDir,'data_*.mat'))
            DataFiles.sort()
            QSubFile=os.path.join(JobDir,('QSubFile.sh'))
            Qfid=open(QSubFile,'write')
            for i in DataFiles:
                # find the data chunk
                DataFileName=os.path.split(i)[-1]
                DataChunk=DataFileName.split('_')[-1].split('.')[0]
                JobFileName=os.path.join(JobDir,('job_%s.sh'%(DataChunk)))
                ParamFileName=os.path.join(JobDir,('param_%s.in'%(DataChunk)))
                # create the parameter files
                # remove the .mat from the end of the filename
                
                fid=open(ParamFileName,'write')
                fid.write('InputData=%s'%i[:-4])
                fid.close()
                # Create the job files
                fid=open(JobFileName,'write')
                fid.write('#!/bin/bash\n')
                fid.write('#PBS -S /bin/bash\n\n')
                fid.write('#Choose MCR directory\n')
                fid.write('MCR=/usr/local/matlab-mcr/v82\n\n')
                fid.write('echo "Running on host: `hostname`"\n')
                fid.write('cd $PBS_O_WORKDIR\n')
                fid.write('echo "Current working directory is `pwd`"\n\n')
                fid.write('echo "Starting run at: `date`"\n\n')
                fid.write('/home/steffejr/deploy/run_RunCompiledBootStrap.sh $MCR %s\n'%ParamFileName)
                fid.write('echo "Job finished at: `date`"\n\n')
                fid.close()
                # For this to work the input to this compiled program has to be a parameter file
                Qfid.write('qsub -l walltime=72:00:00:00,mem=3gb %s\n'%JobFileName)
                
                #Str='qsub -l walltime=72:00:00:00,mem=6gb %s'%JobFileName 
                #os.system(Str)
            Qfid.close()
if __name__=='__main__':
    main()        
#Extract('CogRes')
  
        
        
