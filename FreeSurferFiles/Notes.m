http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FreeSurfer
export SUBJECTS_DIR=/share/studies/CogRes/Subjects/P00004200/S0001/T1
reg-feat2anat --feat ECF_r1_P00004200_S0001.feat/ --subject FreeSurfer 


fslswapdim aseg.auto x z -y aseg.autoROT
