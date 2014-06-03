% Write out basic regression maps using the minimum statistics correction
% for multiple comparisons from the permutation results.
%
% Add the ability to perform FDR correction for multiple comparisons when
% using bootstrapping. This requires the BCa Z and p maps to be saved.
%
% Add the ability to test multiple paths. 
%
% Include test for total effects with the mediation. This will need the
% addition of a third dimension to the DIRECT array.
% 
% Add the ability to test multiple paths in the model.
%
% Add a check for whether the data is "good" or not. This would be an
% implicit maskign sort of approach. At the simplest it should check for
% voxels with NaN and check whether there are "bad" subjects or bad
% "voxels."  Then the user may be provided with an option of removing the
% voxel or subject. If a subject/voxel with NaN values is to be retained
% then there needs to be a procedure for interpolating it's value.
%
% Have permutations saved to a file in small chunks instead of the cluster
% split size. This is valuable for processes that get interrupted. 
% This is now a necessity. The memory demands for multiple permutations is
% HUGE. So it is best to write each permutation to a file after completion.