/*
 * MATLAB Compiler: 4.13 (R2010a)
 * Date: Thu Dec 12 14:19:57 2013
 * Arguments: "-B" "macro_default" "-R" "nodisplay" "-R" "-nojvm" "-R"
 * "-singleCompThread" "-m" "-W" "main" "-T" "link:exe" "-v" "-w" "enable" "-d"
 * "../deploy" "subfnVoxelWiseProcessBatch.m" "subfnProcess.m" "bootSE.m" 
 */

#include "mclmcrrt.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_subfnVoxelWiseProcessBatch_session_key[] = {
    '1', 'C', '6', 'C', 'A', 'E', '3', '7', '8', '6', '6', '8', '7', '5', '0',
    'B', '1', '6', '6', '6', '8', '5', 'D', '0', '0', 'F', 'A', 'A', 'D', '5',
    '1', 'C', '7', '8', 'F', '4', 'B', '0', 'D', '0', '9', 'D', '4', 'B', '4',
    'C', 'C', '9', '9', 'E', '8', 'E', '4', 'B', '6', '3', '3', '3', 'E', '8',
    'A', 'D', '9', '7', '5', '9', '2', '6', '2', 'E', '9', '4', '0', '8', '7',
    'D', 'B', 'C', '7', 'B', 'E', '8', '4', '6', '2', '1', '2', '9', 'D', '0',
    '1', 'C', '6', '8', '9', 'B', '0', 'D', '4', 'B', '0', '2', 'D', '8', 'D',
    'A', 'E', 'F', '8', '5', 'F', '6', '6', '9', 'E', '9', '9', 'B', '0', 'F',
    'A', '1', '5', '7', 'D', 'D', '6', '4', '5', '7', '5', '6', '9', 'E', '3',
    'F', '5', '6', '9', '8', '9', '8', '2', 'E', '4', '4', 'C', '1', '8', '2',
    '8', '2', '3', '5', '4', '3', 'B', '0', '3', '4', 'D', '7', '0', '7', '5',
    '4', '0', '0', '0', '9', 'C', 'E', '7', 'E', 'B', '1', '4', 'C', 'D', 'F',
    'A', '6', '3', 'C', 'B', 'D', 'A', '5', '9', '9', 'C', '3', '9', '9', '7',
    '3', '3', 'F', 'E', 'A', '5', '4', 'C', 'F', 'E', '5', 'C', '3', '6', '4',
    'C', '9', 'F', '1', 'C', '9', '0', 'F', 'F', 'B', 'A', '6', '7', '6', 'E',
    'B', 'F', '2', 'B', '3', '7', 'B', '6', 'E', 'C', '6', 'A', 'B', '4', '4',
    'D', 'E', '0', '4', '4', 'D', 'E', '6', '7', '3', '9', '7', 'A', '4', '4',
    '9', '\0'};

const unsigned char __MCC_subfnVoxelWiseProcessBatch_public_key[] = {
    '3', '0', '8', '1', '9', 'D', '3', '0', '0', 'D', '0', '6', '0', '9', '2',
    'A', '8', '6', '4', '8', '8', '6', 'F', '7', '0', 'D', '0', '1', '0', '1',
    '0', '1', '0', '5', '0', '0', '0', '3', '8', '1', '8', 'B', '0', '0', '3',
    '0', '8', '1', '8', '7', '0', '2', '8', '1', '8', '1', '0', '0', 'C', '4',
    '9', 'C', 'A', 'C', '3', '4', 'E', 'D', '1', '3', 'A', '5', '2', '0', '6',
    '5', '8', 'F', '6', 'F', '8', 'E', '0', '1', '3', '8', 'C', '4', '3', '1',
    '5', 'B', '4', '3', '1', '5', '2', '7', '7', 'E', 'D', '3', 'F', '7', 'D',
    'A', 'E', '5', '3', '0', '9', '9', 'D', 'B', '0', '8', 'E', 'E', '5', '8',
    '9', 'F', '8', '0', '4', 'D', '4', 'B', '9', '8', '1', '3', '2', '6', 'A',
    '5', '2', 'C', 'C', 'E', '4', '3', '8', '2', 'E', '9', 'F', '2', 'B', '4',
    'D', '0', '8', '5', 'E', 'B', '9', '5', '0', 'C', '7', 'A', 'B', '1', '2',
    'E', 'D', 'E', '2', 'D', '4', '1', '2', '9', '7', '8', '2', '0', 'E', '6',
    '3', '7', '7', 'A', '5', 'F', 'E', 'B', '5', '6', '8', '9', 'D', '4', 'E',
    '6', '0', '3', '2', 'F', '6', '0', 'C', '4', '3', '0', '7', '4', 'A', '0',
    '4', 'C', '2', '6', 'A', 'B', '7', '2', 'F', '5', '4', 'B', '5', '1', 'B',
    'B', '4', '6', '0', '5', '7', '8', '7', '8', '5', 'B', '1', '9', '9', '0',
    '1', '4', '3', '1', '4', 'A', '6', '5', 'F', '0', '9', '0', 'B', '6', '1',
    'F', 'C', '2', '0', '1', '6', '9', '4', '5', '3', 'B', '5', '8', 'F', 'C',
    '8', 'B', 'A', '4', '3', 'E', '6', '7', '7', '6', 'E', 'B', '7', 'E', 'C',
    'D', '3', '1', '7', '8', 'B', '5', '6', 'A', 'B', '0', 'F', 'A', '0', '6',
    'D', 'D', '6', '4', '9', '6', '7', 'C', 'B', '1', '4', '9', 'E', '5', '0',
    '2', '0', '1', '1', '1', '\0'};

static const char * MCC_subfnVoxelWiseProcessBatch_matlabpath_data[] = 
  { "subfnVoxelWi/", "$TOOLBOXDEPLOYDIR/", "home/jason/matlab/",
    "share/users/js2746_Jason/Scripts/ProcessModelsNeuroImage/FinalCode/",
    "share/apps/SPM/v8/external/fieldtrip/", "$TOOLBOXMATLABDIR/general/",
    "$TOOLBOXMATLABDIR/ops/", "$TOOLBOXMATLABDIR/lang/",
    "$TOOLBOXMATLABDIR/elmat/", "$TOOLBOXMATLABDIR/randfun/",
    "$TOOLBOXMATLABDIR/elfun/", "$TOOLBOXMATLABDIR/specfun/",
    "$TOOLBOXMATLABDIR/matfun/", "$TOOLBOXMATLABDIR/datafun/",
    "$TOOLBOXMATLABDIR/polyfun/", "$TOOLBOXMATLABDIR/funfun/",
    "$TOOLBOXMATLABDIR/sparfun/", "$TOOLBOXMATLABDIR/scribe/",
    "$TOOLBOXMATLABDIR/graph2d/", "$TOOLBOXMATLABDIR/graph3d/",
    "$TOOLBOXMATLABDIR/specgraph/", "$TOOLBOXMATLABDIR/graphics/",
    "$TOOLBOXMATLABDIR/uitools/", "$TOOLBOXMATLABDIR/strfun/",
    "$TOOLBOXMATLABDIR/imagesci/", "$TOOLBOXMATLABDIR/iofun/",
    "$TOOLBOXMATLABDIR/audiovideo/", "$TOOLBOXMATLABDIR/timefun/",
    "$TOOLBOXMATLABDIR/datatypes/", "$TOOLBOXMATLABDIR/verctrl/",
    "$TOOLBOXMATLABDIR/codetools/", "$TOOLBOXMATLABDIR/helptools/",
    "$TOOLBOXMATLABDIR/demos/", "$TOOLBOXMATLABDIR/timeseries/",
    "$TOOLBOXMATLABDIR/hds/", "$TOOLBOXMATLABDIR/guide/",
    "$TOOLBOXMATLABDIR/plottools/", "toolbox/local/",
    "$TOOLBOXMATLABDIR/datamanager/", "toolbox/compiler/", "toolbox/stats/" };

static const char * MCC_subfnVoxelWiseProcessBatch_classpath_data[] = 
  { "" };

static const char * MCC_subfnVoxelWiseProcessBatch_libpath_data[] = 
  { "" };

static const char * MCC_subfnVoxelWiseProcessBatch_app_opts_data[] = 
  { "" };

static const char * MCC_subfnVoxelWiseProcessBatch_run_opts_data[] = 
  { "nodisplay", "-nojvm", "-singleCompThread" };

static const char * MCC_subfnVoxelWiseProcessBatch_warning_state_data[] = 
  { "off:MATLAB:dispatcher:nameConflict" };


mclComponentData __MCC_subfnVoxelWiseProcessBatch_component_data = { 

  /* Public key data */
  __MCC_subfnVoxelWiseProcessBatch_public_key,

  /* Component name */
  "subfnVoxelWiseProcessBatch",

  /* Component Root */
  "",

  /* Application key data */
  __MCC_subfnVoxelWiseProcessBatch_session_key,

  /* Component's MATLAB Path */
  MCC_subfnVoxelWiseProcessBatch_matlabpath_data,

  /* Number of directories in the MATLAB Path */
  41,

  /* Component's Java class path */
  MCC_subfnVoxelWiseProcessBatch_classpath_data,
  /* Number of directories in the Java class path */
  0,

  /* Component's load library path (for extra shared libraries) */
  MCC_subfnVoxelWiseProcessBatch_libpath_data,
  /* Number of directories in the load library path */
  0,

  /* MCR instance-specific runtime options */
  MCC_subfnVoxelWiseProcessBatch_app_opts_data,
  /* Number of MCR instance-specific runtime options */
  0,

  /* MCR global runtime options */
  MCC_subfnVoxelWiseProcessBatch_run_opts_data,
  /* Number of MCR global runtime options */
  3,
  
  /* Component preferences directory */
  "subfnVoxelWi_C71C43719D4F626190B6A32A138C63A6",

  /* MCR warning status data */
  MCC_subfnVoxelWiseProcessBatch_warning_state_data,
  /* Number of MCR warning status modifiers */
  1,

  /* Path to component - evaluated at runtime */
  NULL

};

#ifdef __cplusplus
}
#endif


