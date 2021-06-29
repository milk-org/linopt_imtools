/**
 * @file    linopt_imtools.h
 * @brief   Function prototypes for linear algebra tools
 *
 * CPU-based lineal algebra tools: decomposition, SVD etc...
 *
 */



#ifndef _LINOPTIMTOOLS_H
#define _LINOPTIMTOOLS_H



#include "compute_SVDdecomp.h"
#include "compute_SVDpseudoInverse.h"
#include "image_construct.h"
#include "image_to_vec.h"
#include "image_fitModes.h"
#include "makeCosRadModes.h"
#include "makeCPAmodes.h"
#include "mask_to_pixtable.h"
#include "vec_to_2Dimage.h"



void __attribute__((constructor)) libinit_linopt_imtools();












/* =============================================================================================== */
/* =============================================================================================== */
/** @name 3. CREATE MODES
 *  create mode basis
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */

imageID linopt_imtools_make1Dpolynomials(
    const char *IDout_name,
    long        NBpts,
    long        MaxOrder,
    float       r0pix
);




///@}






/* =============================================================================================== */
/* =============================================================================================== */
/** @name 4. LINEAR DECOMPOSITION
 *
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */




imageID linopt_imtools_image_construct_stream(
    const char *IDmodes_name,
    const char *IDcoeff_name,
    const char *IDout_name
);





errno_t linopt_compute_1Dfit(
    const char *fnamein,
    long        NBpt,
    long        MaxOrder,
    const char *fnameout,
    int         MODE,
    imageID    *outID
);


errno_t linopt_imtools_image_fitModes(
    const char *ID_name,
    const char *IDmodes_name,
    const char *IDmask_name,
    double      SVDeps,
    const char *IDcoeff_name,
    int         reuse,
    imageID    *outIDcoeff
);



double linopt_imtools_match_slow(
    const char *ID_name,
    const char *IDref_name,
    const char *IDmask_name,
    const char *IDsol_name,
    const char *IDout_name
);


double linopt_imtools_match(
    const char *ID_name,
    const char *IDref_name,
    const char *IDmask_name,
    const char *IDsol_name,
    const char *IDout_name
);

///@}





/* =============================================================================================== */
/* =============================================================================================== */
/** @name 5. OPTIMIZATION
 *
 */
///@{
/* =============================================================================================== */
/* =============================================================================================== */


/**
 * @brief Solve for response matrix given a series of input and output
 *
 *  initial value of RM should be best guess
 *  inmask = 0 over input that are known to produce no response
 */
long linopt_compute_linRM_from_inout(const char *IDinput_name,
                                     const char *IDinmask_name, const char *IDoutput_name, const char *IDRM_name);


///@}






#endif
