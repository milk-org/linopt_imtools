/**
 * @file    linopt_imtools.c
 * @brief   linear optimization tools
 *
 * CPU-based lineal algebra tools: decomposition, SVD etc...
 *
 *
 */



/* ================================================================== */
/* ================================================================== */
/*            MODULE INFO                                             */
/* ================================================================== */
/* ================================================================== */

// module default short name
// all CLI calls to this module functions will be <shortname>.<funcname>
// if set to "", then calls use <funcname>
#define MODULE_SHORTNAME_DEFAULT "lintools"

// Module short description
#define MODULE_DESCRIPTION       "Image linear decomposition and optimization tools"










#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <semaphore.h>
#include <sched.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>


#include <time.h>

#include <fitsio.h>


#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_tools/COREMOD_tools.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "statistic/statistic.h"
#include "info/info.h"
#include "linopt_imtools/linopt_imtools.h"
#include "cudacomp/cudacomp.h"

#include "CommandLineInterface/timeutils.h"

#include "compute_SVDpseudoInverse.h"
#include "image_construct.h"
#include "image_to_vec.h"
#include "image_fitModes.h"
#include "mask_to_pixtable.h"



static long NBPARAM;
static long double C0;
// polynomial coeff (degree = 1)
static long double *polycoeff1 = NULL;
// polynomial coeff (degree = 2)
static long double *polycoeff2 = NULL;
static long dfcnt = 0;








/* ================================================================== */
/* ================================================================== */
/*            INITIALIZE LIBRARY                                      */
/* ================================================================== */
/* ================================================================== */

// Module initialization macro in CLIcore.h
// macro argument defines module name for bindings
//
INIT_MODULE_LIB(linopt_imtools)






errno_t linopt_imtools_vec_to_2DImage_cli()
{
    if(
        CLI_checkarg(1, 4) +
        CLI_checkarg(2, 4) +
        CLI_checkarg(3, 4) +
        CLI_checkarg(4, 3) +
        CLI_checkarg(5, 2) +
        CLI_checkarg(6, 2)
        == 0)
    {
        linopt_imtools_vec_to_2DImage(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.string,
            data.cmdargtoken[4].val.string,
            data.cmdargtoken[5].val.numl,
            data.cmdargtoken[6].val.numl
        );

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}




/* =============================================================================================== */
/* =============================================================================================== */
/*                                                                                                 */
/* 3. CREATE MODES                                                                                 */
/*                                                                                                 */
/* =============================================================================================== */
/* =============================================================================================== */


errno_t linopt_imtools_makeCosRadModes_cli()
{
    if(
        CLI_checkarg(1, 3) +
        CLI_checkarg(2, 2) +
        CLI_checkarg(3, 2) +
        CLI_checkarg(4, 1) +
        CLI_checkarg(5, 1)
        == 0)
    {
        linopt_imtools_makeCosRadModes(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.numl,
            data.cmdargtoken[3].val.numl,
            data.cmdargtoken[4].val.numf,
            data.cmdargtoken[5].val.numf
        );

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}


errno_t linopt_imtools_makeCPAmodes_cli()
{
    if(
        CLI_checkarg(1, 3) +
        CLI_checkarg(2, 2) +
        CLI_checkarg(3, 1) +
        CLI_checkarg(4, 1) +
        CLI_checkarg(5, 1) +
        CLI_checkarg(6, 1)
        == 0)
    {
        linopt_imtools_makeCPAmodes(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.numl,
            data.cmdargtoken[3].val.numf,
            data.cmdargtoken[4].val.numf,
            data.cmdargtoken[5].val.numf,
            data.cmdargtoken[6].val.numf,
            1
        );

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}





/* =============================================================================================== */
/* =============================================================================================== */
/*                                                                                                 */
/* 4. LINEAR DECOMPOSITION                                                                         */
/*                                                                                                 */
/* =============================================================================================== */
/* =============================================================================================== */





errno_t linopt_imtools_image_construct_cli()
{
    if(
        CLI_checkarg(1, 4) +
        CLI_checkarg(2, 4) +
        CLI_checkarg(3, 3)
        == 0)
    {
        linopt_imtools_image_construct(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.string,
            NULL
        );

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}



errno_t linopt_imtools_image_construct_stream_cli()
{
    if(
        CLI_checkarg(1, 4) +
        CLI_checkarg(2, 4) +
        CLI_checkarg(3, 4)
        == 0)
    {
        linopt_imtools_image_construct_stream(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.string
        );

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}



errno_t linopt_compute_SVDdecomp_cli()
{
    if(
        CLI_checkarg(1, 4) +
        CLI_checkarg(2, 3) +
        CLI_checkarg(3, 3)
        == 0)
    {
        linopt_compute_SVDdecomp(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.string
        );

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}


errno_t linopt_compute_SVDpseudoInverse_cli()
{
    if(
        CLI_checkarg(1, 4) +
        CLI_checkarg(2, 3) +
        CLI_checkarg(3, 1) +
        CLI_checkarg(4, 2) +
        CLI_checkarg(5, 3)
        == 0)
    {
        linopt_compute_SVDpseudoInverse(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.numf,
            data.cmdargtoken[4].val.numl,
            data.cmdargtoken[5].val.string
        );

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}


errno_t linopt_compute_1Dfit_cli()
{
    if(
        CLI_checkarg(1, 5) +
        CLI_checkarg(2, 2) +
        CLI_checkarg(3, 2) +
        CLI_checkarg(4, 5) +
        CLI_checkarg(5, 2)
        == 0)
    {
        linopt_compute_1Dfit(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.numl,
            data.cmdargtoken[3].val.numl,
            data.cmdargtoken[4].val.string,
            data.cmdargtoken[5].val.numl,
            NULL
        );

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}




/* =============================================================================================== */
/* =============================================================================================== */
/*                                                                                                 */
/* 5. OPTIMIZATION                                                                                 */
/*                                                                                                 */
/* =============================================================================================== */
/* =============================================================================================== */



errno_t linopt_compute_linRM_from_inout_cli()
{
    if(
        CLI_checkarg(1, 4) +
        CLI_checkarg(2, 4) +
        CLI_checkarg(3, 4) +
        CLI_checkarg(4, 4)
        == 0)
    {
        linopt_compute_linRM_from_inout(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.string,
            data.cmdargtoken[4].val.string
        );

        return CLICMD_SUCCESS;
    }
    else
    {
        return CLICMD_INVALID_ARG;
    }
}












static errno_t init_module_CLI()
{

    /* =============================================================================================== */
    /* =============================================================================================== */
    /*                                                                                                 */
    /* 1. INITIALIZATION                                                                               */
    /*                                                                                                 */
    /* =============================================================================================== */
    /* =============================================================================================== */



    /* =============================================================================================== */
    /* =============================================================================================== */
    /*                                                                                                 */
    /* 2. CONVERSION                                                                                   */
    /*                                                                                                 */
    /* =============================================================================================== */
    /* =============================================================================================== */

    CLIADDCMD_linopt_imtools__mask_to_pixtable();

    CLIADDCMD_linopt_imtools__image_to_vec();


    RegisterCLIcommand(
        "vec2im",
        __FILE__, linopt_imtools_vec_to_2DImage_cli,
        "remap vector to image",
        "<vecname> <pixindex> <pixmult> <imname> <xsize> <ysize>",
        "im2vec vecim pixi pixm im 512 512",
        "long linopt_imtools_vec_to_2DImage(const char *IDvec_name, const char *IDpixindex_name, const char *IDpixmult_name, const char *ID_name, long xsize, long ysize)");



    /* =============================================================================================== */
    /* =============================================================================================== */
    /*                                                                                                 */
    /* 3. CREATE MODES                                                                                 */
    /*                                                                                                 */
    /* =============================================================================================== */
    /* =============================================================================================== */

    RegisterCLIcommand(
        "mkcosrmodes",
        __FILE__,
        linopt_imtools_makeCosRadModes_cli,
        "make basis of cosine radial modes",
        "<output image name> <image size [long]> <kmax [long]> <radius [float]> <overfill factor [float]>",
        "mkcosrmodes cmodes 256 100 80.0 2.0",
        "long linopt_imtools_makeCosRadModes(const char *ID_name, long size, long kmax, float radius, float radfactlim, int writeMfile)");


    RegisterCLIcommand(
        "mkFouriermodes",
        __FILE__,
        linopt_imtools_makeCPAmodes_cli,
        "make basis of Fourier Modes",
        "<output image name> <image size> <CPAmax float> <deltaCPA float> <beam radius> <overfill factor>",
        "mkFouriermodes fmodes 256 10.0 0.8 80.0 2.0",
        "long linopt_imtools_makeCPAmodes(const char *ID_name, long size, float CPAmax, float deltaCPA, float radius, float radfactlim)");



    /* =============================================================================================== */
    /* =============================================================================================== */
    /*                                                                                                 */
    /* 4. LINEAR DECOMPOSITION                                                                         */
    /*                                                                                                 */
    /* =============================================================================================== */
    /* =============================================================================================== */

    CLIADDCMD_linopt_imtools__image_fitModes();

    CLIADDCMD_linopt_imtools__image_construct();




    RegisterCLIcommand(
        "imlinconstructs",
        __FILE__,
        linopt_imtools_image_construct_stream_cli,
        "construct image as linear sum of modes (stream mode)",
        "<modes> <coeffs> <outim>", "imlinconstructs modes coeffs outim",
        "long linopt_imtools_image_construct_stream(const char *IDmodes_name, const char *IDcoeff_name, const char *IDout_name)");


    RegisterCLIcommand(
        "imsvd",
        __FILE__,
        linopt_compute_SVDdecomp_cli,
        "Singular values decomposition",
        "<image cube> <SVD modes> <coeffs>",
        "imsvd imc svdm coeffs",
        "long linopt_compute_SVDdecomp(const char *IDin_name, const char *IDout_name, const char *IDcoeff_name)");


    CLIADDCMD_linopt_imtools__compute_SVDpseudoinverse();


    RegisterCLIcommand(
        "linopt1Dfit",
        __FILE__,
        linopt_compute_1Dfit_cli,
        "least-square 1D fit",
        "<output data file> <NBpt> <fit order> <output coeff file> <fit MODE>",
        "linopt1Dfit data.txt 1000 10 fitsol.txt 0",
        "long linopt_compute_1Dfit(const char *fnamein, long NBpt, long MaxOrder, const char *fnameout, int MODE)");


    /* =============================================================================================== */
    /* =============================================================================================== */
    /*                                                                                                 */
    /* 5. OPTIMIZATION                                                                                 */
    /*                                                                                                 */
    /* =============================================================================================== */
    /* =============================================================================================== */


    RegisterCLIcommand(
        "lincRMiter",
        __FILE__,
        linopt_compute_linRM_from_inout_cli,
        "estimate response matrix from input and output",
        "<input cube> <inmask> <output cube> <RM>",
        "lincRMiter inC inmask outC imRM",
        "long linopt_compute_linRM_iter(const char *IDinput_name, const char *IDinmask_name, const char *IDoutput_name, const char *IDRM_name)");


    // add atexit functions here


    return RETURN_SUCCESS;
}













imageID linopt_imtools_vec_to_2DImage(
    const char *IDvec_name,
    const char *IDpixindex_name,
    const char *IDpixmult_name,
    const char *ID_name,
    long        xsize,
    long        ysize
)
{
    DEBUG_TRACE_FSTART();

    imageID ID;
    imageID IDvec;
    long k;
    imageID IDpixindex, IDpixmult;
    long NBpix;

    IDvec = image_ID(IDvec_name);
    IDpixindex = image_ID(IDpixindex_name);
    IDpixmult = image_ID(IDpixmult_name);
    NBpix = data.image[IDpixindex].md[0].nelement;

    FUNC_CHECK_RETURN(
        create_2Dimage_ID(ID_name, xsize, ysize, &ID));

    for(k = 0; k < NBpix; k++)
    {
        data.image[ID].array.F[data.image[IDpixindex].array.SI64[k]] =
            data.image[IDvec].array.F[k] / data.image[IDpixmult].array.F[k];
    }

    DEBUG_TRACE_FEXIT();
    return ID;
}








/* =============================================================================================== */
/* =============================================================================================== */
/*                                                                                                 */
/* 3. CREATE MODES                                                                                 */
/*                                                                                                 */
/* =============================================================================================== */
/* =============================================================================================== */



// r0pix is r=1 in pixel unit

imageID linopt_imtools_make1Dpolynomials(
    const char *IDout_name,
    long        NBpts,
    long        MaxOrder,
    float       r0pix
)
{
    DEBUG_TRACE_FSTART();

    imageID IDout;
    long xsize, ysize, zsize;
    long ii, kk;

    xsize = NBpts;
    ysize = 1;
    zsize = MaxOrder;

    FUNC_CHECK_RETURN(
        create_3Dimage_ID(IDout_name, xsize, ysize, zsize, &IDout));

    for(kk = 0; kk < zsize; kk++)
    {
        for(ii = 0; ii < xsize; ii++)
        {
            float r = 1.0 * ii / r0pix;
            data.image[IDout].array.F[kk * xsize + ii] = pow(r, 1.0 * kk);
        }
    }

    DEBUG_TRACE_FEXIT();
    return IDout;
}



//
// make cosine radial modes
//
imageID linopt_imtools_makeCosRadModes(
    const char *ID_name,
    long        size,
    long        kmax,
    float       radius,
    float       radfactlim
)
{
    DEBUG_TRACE_FSTART();

    imageID ID;
    long ii, jj;
    long k;
    long size2;
    imageID IDr;
    FILE *fp;

    size2 = size * size;
    create_2Dimage_ID("linopt_tmpr", size, size, &IDr);

    fp = fopen("ModesExpr_CosRad.txt", "w");
    fprintf(fp, "# unit for r = %f pix\n", radius);
    fprintf(fp, "\n");
    for(k = 0; k < kmax; k++)
    {
        fprintf(fp, "%5ld   cos(r*M_PI*%ld)\n", k, k);
    }


    fclose(fp);


    for(ii = 0; ii < size; ii++)
    {
        float x = (1.0 * ii - 0.5 * size) / radius;
        for(jj = 0; jj < size; jj++)
        {
            float y = (1.0 * jj - 0.5 * size) / radius;
            float r = sqrt(x * x + y * y);
            data.image[IDr].array.F[jj * size + ii] = r;
        }
    }

    FUNC_CHECK_RETURN(
        create_3Dimage_ID(ID_name, size, size, kmax, &ID));

    for(k = 0; k < kmax; k++)
        for(ii = 0; ii < size2; ii++)
        {
            float r = data.image[IDr].array.F[ii];
            if(r < radfactlim)
            {
                data.image[ID].array.F[k * size2 + ii] = cos(r * M_PI * k);
            }
        }


    delete_image_ID("linopt_tmpr", DELETE_IMAGE_ERRMODE_WARNING);

    DEBUG_TRACE_FEXIT();

    return ID;
}



long linopt_imtools_makeCPAmodes(
    const char *ID_name,
    long        size,
    float       CPAmax,
    float       deltaCPA,
    float       radius,
    float       radfactlim,
    int         writeMfile
)
{
    imageID ID;
    imageID IDx, IDy, IDr;
    float CPAx, CPAy;
    float x, y, r;
    long ii, jj;
    long k, k1;
    long NBmax;
    float *CPAxarray;
    float *CPAyarray;
    float *CPArarray;
    long size2;
    long NBfrequ;
    //float y0;
    //float ydist;
    float eps;
    FILE *fp;

    long IDfreq;

    eps = 0.1 * deltaCPA;
    printf("size       = %ld\n", size);
    printf("CPAmax     = %f\n", CPAmax);
    printf("deltaCPA   = %f\n", deltaCPA);
    printf("radius     = %f\n", radius);
    printf("radfactlim = %f\n", radfactlim);


    size2 = size * size;
    create_2Dimage_ID("cpa_tmpx", size, size, &IDx);
    create_2Dimage_ID("cpa_tmpy", size, size, &IDy);
    create_2Dimage_ID("cpa_tmpr", size, size, &IDr);


    printf("precomputing x, y, r\n");
    fflush(stdout);

    for(ii = 0; ii < size; ii++)
    {
        x = (1.0 * ii - 0.5 * size) / radius;
        for(jj = 0; jj < size; jj++)
        {
            y = (1.0 * jj - 0.5 * size) / radius;
            r = sqrt(x * x + y * y);
            data.image[IDx].array.F[jj * size + ii] = x;
            data.image[IDy].array.F[jj * size + ii] = y;
            data.image[IDr].array.F[jj * size + ii] = r;
        }
    }


    printf("CPA: max = %f   delta = %f\n", CPAmax, deltaCPA);
    fflush(stdout);
    NBfrequ = 0;
    for(CPAx = 0; CPAx < CPAmax; CPAx += deltaCPA)
        for(CPAy = -CPAmax; CPAy < CPAmax; CPAy += deltaCPA)
        {
            NBfrequ ++;
        }

    printf("NBfrequ = %ld\n", NBfrequ);
    fflush(stdout);

    CPAxarray = (float *) malloc(sizeof(float) * NBfrequ);
    if(CPAxarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    CPAyarray = (float *) malloc(sizeof(float) * NBfrequ);
    if(CPAyarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    CPArarray = (float *) malloc(sizeof(float) * NBfrequ);
    if(CPArarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    NBfrequ = 0;
    //ydist = 2.0*deltaCPA;
    //y0 = 0.0;
    for(CPAx = 0; CPAx < CPAmax; CPAx += deltaCPA)
    {
        for(CPAy = 0; CPAy < CPAmax; CPAy += deltaCPA)
        {
            CPAxarray[NBfrequ] = CPAx;
            CPAyarray[NBfrequ] = CPAy;
            CPArarray[NBfrequ] = sqrt(CPAx * CPAx + CPAy * CPAy);
            NBfrequ++;
        }
        if(CPAx > eps)
        {
            for(CPAy = -deltaCPA; CPAy > -CPAmax; CPAy -= deltaCPA)
            {
                CPAxarray[NBfrequ] = CPAx;
                CPAyarray[NBfrequ] = CPAy;
                CPArarray[NBfrequ] = sqrt(CPAx * CPAx + CPAy * CPAy);
                NBfrequ++;
            }
        }
    }


    //  for(k1=0;k1<NBfrequ;k1++)
    //printf("%ld %f %f %f\n", k1, CPAxarray[k1], CPAyarray[k1], CPArarray[k1]);



    //  printf("sorting\n");
    // fflush(stdout);

    quick_sort3_float(CPArarray, CPAxarray, CPAyarray, NBfrequ);




    NBmax = NBfrequ * 2;
    create_3Dimage_ID(ID_name, size, size, NBmax - 1, &ID);



    if(writeMfile == 1)
    {
        fp = fopen("ModesExpr_CPA.txt", "w");
        fprintf(fp, "# size       = %ld\n", size);
        fprintf(fp, "# CPAmax     = %f\n", CPAmax);
        fprintf(fp, "# deltaCPA   = %f\n", deltaCPA);
        fprintf(fp, "# radius     = %f\n", radius);
        fprintf(fp, "# radfactlim = %f\n", radfactlim);
        fprintf(fp, "# \n");
        fprintf(fp, "# Unit for x and y = radius [pixel]\n");
        fprintf(fp, "# \n");
        fprintf(fp, "%4ld %10.5f %10.5f    1.0\n", (long) 0, 0.0, 0.0);
        k1 = 1;
        k = 2;
        while(k < NBmax)
        {
            CPAx = CPAxarray[k1];
            CPAy = CPAyarray[k1];
            if(CPAy < 0)
            {
                fprintf(fp, "%4ld %10.5f %10.5f    cos(M_PI*(x*%.5f-y*%.5f))\n", k - 1, CPAx,
                        CPAy, CPAx, -CPAy);
                fprintf(fp, "%4ld %10.5f %10.5f    sin(M_PI*(x*%.5f-y*%.5f))\n", k, CPAx, CPAy,
                        CPAx, -CPAy);
            }
            else
            {
                fprintf(fp, "%4ld %10.5f %10.5f    cos(M_PI*(x*%.5f+y*%.5f))\n", k - 1, CPAx,
                        CPAy, CPAx, CPAy);
                fprintf(fp, "%4ld %10.5f %10.5f    sin(M_PI*(x*%.5f+y*%.5f))\n", k, CPAx, CPAy,
                        CPAx, CPAy);
            }
            k += 2;
            k1++;
        }

        fclose(fp);
    }

    delete_image_ID("cpamodesfreq", DELETE_IMAGE_ERRMODE_WARNING);
    create_2Dimage_ID("cpamodesfreq", NBmax - 1, 1, &IDfreq);


    // mode 0 (piston)
    data.image[IDfreq].array.F[0] = 0.0;
    for(ii = 0; ii < size2; ii++)
    {
        x = data.image[IDx].array.F[ii];
        y = data.image[IDy].array.F[ii];
        r = data.image[IDr].array.F[ii];
        if(r < radfactlim)
        {
            data.image[ID].array.F[ii] = 1.0;
        }
    }

    k1 = 1;
    k = 2;
    while(k < NBmax)
    {
        //      printf("\r%5ld / %5ld          ", k, NBmax);
        //      fflush(stdout);
        CPAx = CPAxarray[k1];
        CPAy = CPAyarray[k1];
        // printf("    %ld %f %f\n", k1, CPAx, CPAy);
        for(ii = 0; ii < size2; ii++)
        {
            x = data.image[IDx].array.F[ii];
            y = data.image[IDy].array.F[ii];
            r = data.image[IDr].array.F[ii];
            data.image[IDfreq].array.F[k - 1] = sqrt(CPAx * CPAx + CPAy * CPAy);
            data.image[IDfreq].array.F[k] = sqrt(CPAx * CPAx + CPAy * CPAy);
            if(r < radfactlim)
            {
                data.image[ID].array.F[(k - 1)*size2 + ii] = cos(M_PI * (x * CPAx + y * CPAy));
                data.image[ID].array.F[k * size2 + ii] = sin(M_PI * (x * CPAx + y * CPAy));
            }
        }
        k += 2;
        k1++;
    }
    //  printf("done \n");
    // fflush(stdout);

    free(CPAxarray);
    free(CPAyarray);
    free(CPArarray);


    delete_image_ID("cpa_tmpx",DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("cpa_tmpy", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("cpa_tmpr", DELETE_IMAGE_ERRMODE_WARNING);

    // printf("done \n");
    //fflush(stdout);



    return NBmax;
}






/* =============================================================================================== */
/* =============================================================================================== */
/*                                                                                                 */
/* 4. LINEAR DECOMPOSITION                                                                         */
/*                                                                                                 */
/* =============================================================================================== */
/* =============================================================================================== */





/* --------------------------------------------------------------- */
/*                                                                 */
/*           Functions for optimization                            */
/*                                                                 */
/* --------------------------------------------------------------- */

double linopt_imtools_opt_f(
    const gsl_vector *v,
    __attribute__((unused)) void *params
)
{
    double value;
    long k, l, n;


    n = NBPARAM;
    value = C0;
    for(k = 0; k < n; k++)
    {
        value += polycoeff1[k] * gsl_vector_get(v, k);
    }
    for(k = 0; k < n; k++)
        for(l = 0; l < n; l++)
        {
            value += polycoeff2[l * n + k] * gsl_vector_get(v, k) * gsl_vector_get(v, l);
        }

    return(value);
}




void linopt_imtools_opt_df(
    const gsl_vector *v,
    void             *params,
    gsl_vector       *df
)
{
    double epsilon = 1.0e-8;
    long i, j;
    double v1, v2;
    gsl_vector *vcp;

    vcp = gsl_vector_alloc(NBPARAM);
    //v0 = linopt_imtools_opt_f (v, params);

    for(i = 0; i < NBPARAM; i++)
    {
        for(j = 0; j < NBPARAM; j++)
        {
            gsl_vector_set(vcp, j, gsl_vector_get(v, j));
        }
        gsl_vector_set(vcp, i, gsl_vector_get(v, i) + epsilon);
        v1 = linopt_imtools_opt_f(vcp, params);
        gsl_vector_set(vcp, i, gsl_vector_get(v, i) - epsilon);
        v2 = linopt_imtools_opt_f(vcp, params);
        gsl_vector_set(df, i, (double)((v1 - v2) / (2.0 * epsilon)));
    }

    if(0)
    {
        printf("%ld df = (", dfcnt);
        for(i = 0; i < NBPARAM; i++)
        {
            printf(" %g", gsl_vector_get(df, i));
        }
        printf(" )\n");
    }
    dfcnt ++;

    if(dfcnt > 50)
    {
        exit(0);
    }

    gsl_vector_free(vcp);
}






void linopt_imtools_opt_fdf(
    const gsl_vector *x,
    void             *params,
    double           *f,
    gsl_vector       *df
)
{
    *f = linopt_imtools_opt_f(x, params);
    linopt_imtools_opt_df(x, params, df);
}










// FLOAT only
imageID linopt_imtools_image_construct_stream(
    const char *IDmodes_name,
    const char *IDcoeff_name,
    const char *IDout_name
)
{
    imageID IDout;
    imageID IDmodes;
    imageID IDcoeff;
    long ii, kk;
    long xsize, ysize, zsize;
    long sizexy;
    int semval;
    uint64_t cnt = 0;
    int RT_priority = 80; //any number from 0-99
    struct sched_param schedpar;
    int NOSEM = 1; // ignore input semaphore, use counter


    schedpar.sched_priority = RT_priority;
#ifndef __MACH__
    sched_setscheduler(0, SCHED_FIFO,
                       &schedpar); //other option is SCHED_RR, might be faster
#endif



    IDmodes = image_ID(IDmodes_name);
    //datatype = data.image[IDmodes].md[0].datatype;

    xsize = data.image[IDmodes].md[0].size[0];
    ysize = data.image[IDmodes].md[0].size[1];
    zsize = data.image[IDmodes].md[0].size[2];

    sizexy = xsize * ysize;

    if(variable_ID("NOSEM") != -1)
    {
        NOSEM = 1;
    }
    else
    {
        NOSEM = 0;
    }

    IDout = image_ID(IDout_name);
    IDcoeff = image_ID(IDcoeff_name);

    while(1 == 1)
    {
        if((data.image[IDcoeff].md[0].sem == 0) || (NOSEM == 1))
        {
            while(cnt == data.image[IDcoeff].md[0].cnt0) // test if new frame exists
            {
                usleep(5);
            }
            cnt = data.image[IDcoeff].md[0].cnt0;
        }
        else
        {
            sem_wait(data.image[IDcoeff].semptr[0]);
        }

        for(ii = 0; ii < sizexy; ii++)
        {
            data.image[IDout].array.F[ii] = 0.0;
        }

        data.image[IDout].md[0].write = 1;
        for(kk = 0; kk < zsize; kk++)
            for(ii = 0; ii < sizexy; ii++)
            {
                data.image[IDout].array.F[ii] += data.image[IDcoeff].array.F[kk] *
                                                 data.image[IDmodes].array.F[kk * sizexy + ii];
            }
        sem_getvalue(data.image[IDout].semptr[0], &semval);
        if(semval < SEMAPHORE_MAXVAL)
        {
            sem_post(data.image[IDout].semptr[0]);
        }

        data.image[IDout].md[0].cnt0++;
        data.image[IDout].md[0].write = 0;
    }

    return IDout;
}




// rotation matrix written as SVD_VTm

imageID linopt_compute_SVDdecomp(
    const char *IDin_name,
    const char *IDout_name,
    const char *IDcoeff_name
)
{
    DEBUG_TRACE_FSTART();

    imageID IDin;
    imageID IDout;
    imageID IDcoeff;
    long k, ii;
    gsl_matrix *matrix_D; /* input */
    gsl_matrix *matrix_Dtra;
    gsl_matrix *matrix_DtraD;
    gsl_matrix *matrix_DtraD_evec;
    //gsl_matrix *matrix1;
    //gsl_matrix *matrix2;
    gsl_vector *matrix_DtraD_eval;
    gsl_eigen_symmv_workspace *w;
    gsl_matrix *matrix_save;


    long m;
    long n;
    uint32_t *arraysizetmp;

    long kk, kk1;
    imageID ID_VTmatrix;

    arraysizetmp = (uint32_t *) malloc(sizeof(uint32_t) * 3);
    if(arraysizetmp == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }


    printf("[SVD start]");
    fflush(stdout);


    IDin = image_ID(IDin_name);


    n = data.image[IDin].md[0].size[0] * data.image[IDin].md[0].size[1];
    m = data.image[IDin].md[0].size[2];

    matrix_DtraD_eval = gsl_vector_alloc(m);
    matrix_D = gsl_matrix_alloc(n, m);
    matrix_Dtra = gsl_matrix_alloc(m, n);
    matrix_DtraD = gsl_matrix_alloc(m, m);
    matrix_DtraD_evec = gsl_matrix_alloc(m, m);


    /* write matrix_D */
    for(k = 0; k < m; k++)
    {
        for(ii = 0; ii < n; ii++)
        {
            gsl_matrix_set(matrix_D, ii, k, data.image[IDin].array.F[k * n + ii]);
        }
    }
    /* compute DtraD */
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, matrix_D, matrix_D, 0.0,
                   matrix_DtraD);

    /* compute the inverse of DtraD */

    /* first, compute the eigenvalues and eigenvectors */
    w =   gsl_eigen_symmv_alloc(m);
    matrix_save = gsl_matrix_alloc(m, m);
    gsl_matrix_memcpy(matrix_save, matrix_DtraD);
    gsl_eigen_symmv(matrix_save, matrix_DtraD_eval, matrix_DtraD_evec, w);

    gsl_matrix_free(matrix_save);
    gsl_eigen_symmv_free(w);
    gsl_eigen_symmv_sort(matrix_DtraD_eval, matrix_DtraD_evec,
                         GSL_EIGEN_SORT_ABS_DESC);

    create_2Dimage_ID(IDcoeff_name, m, 1, &IDcoeff);


    for(k = 0; k < m; k++)
    {
        data.image[IDcoeff].array.F[k] = gsl_vector_get(matrix_DtraD_eval, k);
    }



    /** Write rotation matrix to go from DM modes to eigenmodes */
    arraysizetmp[0] = m;
    arraysizetmp[1] = m;
    ID_VTmatrix = image_ID("SVD_VTm");
    if(ID_VTmatrix != -1)
    {
        delete_image_ID("SVD_VTm", DELETE_IMAGE_ERRMODE_WARNING);
    }
    create_image_ID("SVD_VTm", 2, arraysizetmp, _DATATYPE_FLOAT, 0,
                    0, 0, &ID_VTmatrix);
    for(ii = 0; ii < m; ii++) // modes
        for(k = 0; k < m; k++) // modes
        {
            data.image[ID_VTmatrix].array.F[k * m + ii] = (float) gsl_matrix_get(
                        matrix_DtraD_evec, k, ii);
        }

    /// Compute SVD decomp

    create_3Dimage_ID(IDout_name,
                      data.image[IDin].md[0].size[0],
                      data.image[IDin].md[0].size[1],
                      data.image[IDin].md[0].size[2],
                      &IDout);

    for(kk = 0; kk < m; kk++) /// eigen mode index
    {
        //        printf("eigenmode %4ld / %4ld  %g\n", kk, m, data.image[IDcoeff].array.F[kk]);
        //       fflush(stdout);
        for(kk1 = 0; kk1 < m; kk1++)
        {
            for(ii = 0; ii < n; ii++)
            {
                data.image[IDout].array.F[kk * n + ii] += data.image[ID_VTmatrix].array.F[kk1 *
                        m + kk] * data.image[IDin].array.F[kk1 * n + ii];
            }
        }
    }

    //   delete_image_ID("SVD_VTm");


    free(arraysizetmp);

    gsl_matrix_free(matrix_D);
    gsl_matrix_free(matrix_Dtra);
    gsl_matrix_free(matrix_DtraD);
    gsl_matrix_free(matrix_DtraD_evec);
    gsl_vector_free(matrix_DtraD_eval);

    printf("[SVD done]\n");
    fflush(stdout);

    DEBUG_TRACE_FEXIT();
    return IDout;
}






// MODE :
// 0 : polynomial
//
errno_t linopt_compute_1Dfit(
    const char *fnamein,
    long        NBpt,
    long        MaxOrder,
    const char *fnameout,
    int         MODE,
    imageID    *outID
)
{
    DEBUG_TRACE_FSTART();

    float *xarray;
    float *valarray;

    FILE *fp;
    long ii;

    imageID IDin, IDin0;
    imageID IDmask;
    imageID IDmodes;
    long NBmodes;
    long m;

    float SVDeps = 0.0000001;

    long IDout, IDout0;
    double val, vale, err;

    long NBiter = 100;
    float gain = 1.0;
    long iter;


    xarray = (float *) malloc(sizeof(float) * NBpt);
    if(xarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    valarray = (float *) malloc(sizeof(float) * NBpt);
    if(valarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }


    fp = fopen(fnamein, "r");
    for(ii = 0; ii < NBpt; ii++)
    {
        int fscanfcnt = fscanf(fp, "%f %f\n", &xarray[ii], &valarray[ii]);

        if(fscanfcnt == EOF)
        {
            if(ferror(fp))
            {
                perror("fscanf");
            }
            else
            {
                fprintf(stderr,
                        "Error: fscanf reached end of file, no matching characters, no matching failure\n");
            }
            exit(EXIT_FAILURE);
        }
        else if(fscanfcnt != 2)
        {
            fprintf(stderr,
                    "Error: fscanf successfully matched and assigned %i input items, 2 expected\n",
                    fscanfcnt);
            exit(EXIT_FAILURE);
        }
    }
    fclose(fp);

    create_2Dimage_ID("invect", NBpt, 1, &IDin);
    create_2Dimage_ID("invect0", NBpt, 1, &IDin0);
    create_2Dimage_ID("inmask", NBpt, 1, &IDmask);

    for(ii = 0; ii < NBpt; ii++)
    {
        //			printf("%18.16f  %+18.16f\n", xarray[ii], valarray[ii]);
        data.image[IDin].array.F[ii] = valarray[ii];
        data.image[IDin0].array.F[ii] = valarray[ii];
        data.image[IDmask].array.F[ii] = 1.0;
    }

    NBmodes = MaxOrder;
    create_3Dimage_ID("fitmodes", NBpt, 1, NBmodes, &IDmodes);
    create_2Dimage_ID("outcoeff", NBmodes, 1, &IDout);

    switch(MODE)
    {
    case 0 :
        for(m = 0; m < NBmodes; m++)
        {
            for(ii = 0; ii < NBpt; ii++)
            {
                data.image[IDmodes].array.F[m * NBpt + ii] = pow(xarray[ii], 1.0 * m);
            }
        }
        break;
    case 1 :
        for(m = 0; m < NBmodes; m++)
        {
            for(ii = 0; ii < NBpt; ii++)
            {
                data.image[IDmodes].array.F[m * NBpt + ii] = cos(xarray[ii] * M_PI * m);
            }
        }
        break;
    default :
        printf("ERROR: MODE = %d not supported\n", MODE);
        exit(0);
        break;
    }

    list_image_ID();

    for(iter = 0; iter < NBiter; iter++)
    {
        FUNC_CHECK_RETURN(
            linopt_imtools_image_fitModes("invect0", "fitmodes", "inmask", SVDeps,
                                          "outcoeffim0", 1, NULL));
        IDout0 = image_ID("outcoeffim0");


        for(m = 0; m < NBmodes; m++)
        {
            data.image[IDout].array.F[m] += gain * data.image[IDout0].array.F[m];
        }

        for(ii = 0; ii < NBpt; ii++)
        {
            err = 0.0;
            val = 0.0;
            for(m = 0; m < NBmodes; m++)
            {
                val += data.image[IDout].array.F[m] * data.image[IDmodes].array.F[m * NBpt +
                        ii];
            }
            data.image[IDin0].array.F[ii] = data.image[IDin].array.F[ii] - val;
            err += data.image[IDin0].array.F[ii] * data.image[IDin0].array.F[ii];
        }
        err = sqrt(err / NBpt);
        printf("ITERATION %4ld   residual = %20g   [gain = %20g]\n", iter, err, gain);
        gain *= 0.95;
    }




    fp = fopen(fnameout, "w");
    for(m = 0; m < NBmodes; m++)
    {
        fprintf(fp, "%4ld %+.8g\n", m, data.image[IDout].array.F[m]);
    }
    fclose(fp);


    fp = fopen("testout.txt", "w");
    err = 0.0;
    for(ii = 0; ii < NBpt; ii++)
    {
        val = 0.0;
        for(m = 0; m < NBmodes; m++)
        {
            val += data.image[IDout].array.F[m] * data.image[IDmodes].array.F[m * NBpt +
                    ii];
        }
        vale = valarray[ii] - val;
        err += vale * vale;
        fprintf(fp, "%05ld  %18.16f  %18.16f   %18.16f\n", ii, xarray[ii], valarray[ii],
                val);
    }
    fclose(fp);
    err = sqrt(err / NBpt);

    printf("FIT error = %g m\n", err);

    free(xarray);
    free(valarray);

    if(outID != NULL)
    {
        *outID = IDout;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}








//
// match a single image (ID_name) to a linear sum of images within IDref_name
// result is a 1D array of coefficients in IDsol_name
//
double linopt_imtools_match_slow(
    const char *ID_name,
    const char *IDref_name,
    const char *IDmask_name,
    const char *IDsol_name,
    const char *IDout_name
)
{
    long ID, IDref, IDmask, IDsol, IDout;
    long naxes[2];
    long n; // number of reference frames
    long ii, k, l;

    long double val;
    long double valbest;

    // initial random search
    long riter;
    long riterMax = 1000000;

    long double v0;
    long double *tarray = NULL; // temporary array to store values for fixed pixel


    // ref image coefficients (solutions)
    long double *alpha = NULL;
    long double *alphabest = NULL;
    long double ampl;

    /*
      the optimization problem is first rewritten as a 2nd degree polynomial of alpha values
      val = V0 + SUM_{k=0...n-1}{polycoeff1[k]*alpha[k] + SUM_{k=0...n-1}{l=0...k}{polycoeff2[k,l]*alpha[k]*alpha[l]}
     */

    long iter = 0;
    double *params;
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *sminimizer;
    long i;
    gsl_vector *x;
    gsl_multimin_function_fdf opt_func;
    int status;

    //  printf("Input params : %s %s %s\n",ID_name,IDref_name,IDsol_name);

    params = (double *) malloc(sizeof(double) * 1);
    if(params == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    params[0] = 0.0;


    ID = image_ID(ID_name);
    naxes[0] = data.image[ID].md[0].size[0];
    naxes[1] = data.image[ID].md[0].size[1];

    IDmask = image_ID(IDmask_name);
    IDref = image_ID(IDref_name);
    n = data.image[IDref].md[0].size[2];

    printf("Number of points = %ld x %ld\n", naxes[0]*naxes[1], n);


    alpha = (long double *) malloc(sizeof(long double) * n);
    if(alpha == NULL)
    {
        PRINT_ERROR("Cannot allocate memory");
        exit(0);
    }
    alphabest = (long double *) malloc(sizeof(long double) * n);
    if(alphabest == NULL)
    {
        PRINT_ERROR("Cannot allocate memory");
        exit(0);
    }



    polycoeff1 = (long double *) malloc(sizeof(long double) * n);
    if(polycoeff1 == NULL)
    {
        PRINT_ERROR("Cannot allocate memory");
        exit(0);
    }
    polycoeff2 = (long double *) malloc(sizeof(long double) * n * n);
    if(polycoeff2 == NULL)
    {
        PRINT_ERROR("Cannot allocate memory");
        exit(0);
    }

    tarray = (long double *) malloc(sizeof(long double) * n);
    if(tarray == NULL)
    {
        PRINT_ERROR("Cannot allocate memory");
        exit(0);
    }



    // initialize all coeffs to zero
    C0 = 0.0;
    for(k = 0; k < n; k++)
    {
        alpha[k] = 1.0 / n;
        polycoeff1[k] = 0.0;
        for(l = 0; l < n; l++)
        {
            polycoeff2[l * n + k] = 0.0;
        }
    }

    // compute polynomial coefficients
    for(ii = 0; ii < naxes[0]*naxes[1]; ii++)
    {
        v0 = (long double)(data.image[ID].array.F[ii] * data.image[IDmask].array.F[ii]);
        for(k = 0; k < n; k++)
        {
            tarray[k] = (long double)(data.image[IDref].array.F[naxes[0] * naxes[1] * k +
                                      ii] * data.image[IDmask].array.F[ii]);
        }
        C0 += v0 * v0;
        for(k = 0; k < n; k++)
        {
            polycoeff1[k] += -2.0 * v0 * tarray[k];
        }
        for(k = 0; k < n; k++)
            for(l = 0; l < n; l++)
            {
                polycoeff2[l * n + k] += tarray[k] * tarray[l];
            }
    }

    // find solution
    /*   val = C0 + SUM_{k=0...n-1}{polycoeff1[k]*alpha[k] + SUM_{k=0...n-1}{l=0...k}{polycoeff2[k,l]*alpha[k]*alpha[l]}
     */
    val = C0;
    for(k = 0; k < n; k++)
    {
        val += polycoeff1[k] * alpha[k];
    }
    for(k = 0; k < n; k++)
        for(l = 0; l < n; l++)
        {
            val += polycoeff2[l * n + k] * alpha[k] * alpha[l];
        }


    for(k = 0; k < n; k++)
    {
        printf("%g ", (double) alpha[k]);
    }
    printf("-> %g\n", (double) val);
    for(k = 0; k < n; k++)
    {
        alphabest[k] = alpha[k];
    }
    valbest = val;





    for(riter = 0; riter < riterMax; riter++)
    {
        ampl = pow(ran1(), 4.0);
        for(k = 0; k < n; k++)
        {
            alpha[k] = alphabest[k] + ampl * (1.0 - 2.0 * ran1()) / n;
        }

        val = C0;
        for(k = 0; k < n; k++)
        {
            val += polycoeff1[k] * alpha[k];
        }
        for(k = 0; k < n; k++)
            for(l = 0; l < n; l++)
            {
                val += polycoeff2[l * n + k] * alpha[k] * alpha[l];
            }
        if(val < valbest)
        {
            //printf("[%ld/%ld] ",riter,riterMax);
            //for(k=0;k<n;k++)
            //  printf(" %g ", (double) alpha[k]);
            //printf("-> %g\n", (double) val);
            for(k = 0; k < n; k++)
            {
                alphabest[k] = alpha[k];
            }
            valbest = val;
        }
    }

    NBPARAM = n;

    x = gsl_vector_alloc(n);

    for(i = 0; i < n; i++)
    {
        gsl_vector_set(x, i, alphabest[i]);
    }
    printf("Value = %g\n", linopt_imtools_opt_f(x, params));



    opt_func.n = n;
    opt_func.f = &linopt_imtools_opt_f;
    opt_func.df = &linopt_imtools_opt_df;
    opt_func.fdf = &linopt_imtools_opt_fdf;
    opt_func.params = &params;

    x = gsl_vector_alloc(n);

    for(i = 0; i < n; i++)
    {
        gsl_vector_set(x, i, alphabest[i]);
    }

    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    sminimizer = gsl_multimin_fdfminimizer_alloc(T, n);

    gsl_multimin_fdfminimizer_set(sminimizer, &opt_func, x, 1.0e-5, 0.1);

    do
    {
        iter++;
        dfcnt = 0;
        status = gsl_multimin_fdfminimizer_iterate(sminimizer);
        if(status)
        {
            break;
        }
        status = gsl_multimin_test_gradient(sminimizer->gradient, 1e-5);
        if(status == GSL_SUCCESS)
        {
            printf("Minimum found at:\n");
            printf("%5ld : ", iter);
            //for(i=0;i<n;i++)
            // printf("%.8f ",gsl_vector_get(sminimizer->x, i));
            printf("    %10.8f\n", sminimizer->f);
        }
    }
    while(status == GSL_CONTINUE && iter < 1000);

    for(i = 0; i < n; i++)
    {
        alphabest[i] = gsl_vector_get(sminimizer->x, i);
    }

    for(i = 0; i < n; i++)
    {
        gsl_vector_set(x, i, alphabest[i]);
    }
    printf("Value after minimization = %g\n", linopt_imtools_opt_f(x, params));

    gsl_multimin_fdfminimizer_free(sminimizer);
    gsl_vector_free(x);


    create_2Dimage_ID(IDsol_name, n, 1, &IDsol);
    for(i = 0; i < n; i++)
    {
        data.image[IDsol].array.F[i] = alphabest[i];
    }



    // compute residual

    create_2Dimage_ID(IDout_name, naxes[0], naxes[1], &IDout);

    for(ii = 0; ii < naxes[0]*naxes[1]; ii++)
    {
        data.image[IDout].array.F[ii] = 0.0;
    }
    for(k = 0; k < n; k++)
        for(ii = 0; ii < naxes[0]*naxes[1]; ii++)
        {
            data.image[IDout].array.F[ii] += alphabest[k] *
                                             data.image[IDref].array.F[naxes[0] * naxes[1] * k + ii];
        }


    free(alpha);
    alpha = NULL;
    free(alphabest);
    alphabest = NULL;
    free(polycoeff1);
    polycoeff1 = NULL;
    free(polycoeff2);
    polycoeff2 = NULL;
    free(tarray);
    tarray = NULL;

    free(params);

    return((double) val);
}





// match a single image (ID_name) to a linear sum of images within IDref_name
// result is a 1D array of coefficients in IDsol_name
//
// n = number of observations
// p = number of variables
//
// ID_name is input, size (n,1)
// IDsol_name must contain initial solution
//

double linopt_imtools_match(
    const char *ID_name,
    const char *IDref_name,
    const char *IDmask_name,
    const char *IDsol_name,
    const char *IDout_name
)
{
    gsl_multifit_linear_workspace *work;
    uint32_t n, p;
    imageID ID, IDref, IDmask, IDsol, IDout;
    long naxes[3];
    long i, j, k, ii;
    gsl_matrix *X;
    gsl_vector *y; // measurements
    gsl_vector *c;
    gsl_vector *w;
    gsl_matrix *cov;
    double chisq;

    ID = image_ID(ID_name);
    naxes[0] = data.image[ID].md[0].size[0];
    naxes[1] = data.image[ID].md[0].size[1];
    n = naxes[0] * naxes[1];
    IDmask = image_ID(IDmask_name);
    IDref = image_ID(IDref_name);
    p = data.image[IDref].md[0].size[2];

    // some verification
    if(IDref == -1)
    {
        PRINT_ERROR("input ref missing\n");
    }
    if(IDmask == -1)
    {
        PRINT_ERROR("input mask missing\n");
    }
    if(data.image[IDmask].md[0].size[0] != data.image[ID].md[0].size[0])
    {
        PRINT_ERROR("mask size[0] is wrong\n");
    }
    if(data.image[IDmask].md[0].size[1] != data.image[ID].md[0].size[1])
    {
        PRINT_ERROR("mask size[1] is wrong\n");
    }


    printf("n,p = %ld %ld\n", (long) n, (long) p);
    fflush(stdout);

    y = gsl_vector_alloc(n);  // measurements
    for(i = 0; i < n; i++)
    {
        gsl_vector_set(y, i, data.image[ID].array.F[i]);
    }

    w = gsl_vector_alloc(n);
    for(i = 0; i < n; i++)
    {
        gsl_vector_set(w, i, data.image[IDmask].array.F[i]);
    }

    X = gsl_matrix_alloc(n, p);
    for(i = 0; i < n; i++)
        for(j = 0; j < p; j++)
        {
            gsl_matrix_set(X, i, j, data.image[IDref].array.F[j * n + i]);
        }
    c = gsl_vector_alloc(p);  // solution (coefficients)
    cov = gsl_matrix_alloc(p, p);

    work = gsl_multifit_linear_alloc(n, p);

    printf("-");
    fflush(stdout);
    gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
    printf("-");
    fflush(stdout);

    create_2Dimage_ID(IDsol_name, p, 1, &IDsol);
    for(i = 0; i < p; i++)
    {
        data.image[IDsol].array.F[i] = gsl_vector_get(c, i);
    }

    gsl_multifit_linear_free(work);
    gsl_vector_free(y);
    gsl_vector_free(w);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_matrix_free(cov);

    printf(" . ");
    fflush(stdout);

    // compute residual
    create_2Dimage_ID(IDout_name, naxes[0], naxes[1], &IDout);
    for(ii = 0; ii < naxes[0]*naxes[1]; ii++)
    {
        data.image[IDout].array.F[ii] = 0.0;
    }
    for(k = 0; k < p; k++)
        for(ii = 0; ii < naxes[0]*naxes[1]; ii++)
        {
            data.image[IDout].array.F[ii] += data.image[IDsol].array.F[k] *
                                             data.image[IDref].array.F[naxes[0] * naxes[1] * k + ii];
        }

    return(chisq);
}













/* =============================================================================================== */
/* =============================================================================================== */
/*                                                                                                 */
/* 5. OPTIMIZATION                                                                                 */
/*                                                                                                 */
/* =============================================================================================== */
/* =============================================================================================== */



//
// solve for response matrix given a series of input and output
// initial value of RM should be best guess
// inmask = 0 over input that are known to produce no response
//
imageID linopt_compute_linRM_from_inout(
    const char *IDinput_name,
    const char *IDinmask_name,
    const char *IDoutput_name,
    const char *IDRM_name
)
{
    imageID IDRM;
    imageID IDin;
    imageID IDinmask;
    imageID IDout;
    long insize; // number of input
    long xsizein, ysizein, xsizeout, ysizeout;
    double fitval;
    long kk, ii_in, jj_in, ii_out, jj_out;
    //double tot;
    imageID IDtmp;
    double tmpv1;
    //long iter;
    imageID IDout1;
    //double alpha = 0.001;

    uint32_t *sizearray;
    imageID IDpokeM; // poke matrix (input)
    //imageID IDoutM; // outputX
    double SVDeps = 1.0e-4;

    long NBact, act;
    long *inpixarray;
    long spl; // sample measurement
    long ii;
    imageID ID_rm;
    int autoMask_MODE =
        0; // if 1, automatically measure input mask based on IDinput_name image
    imageID IDpinv;
    //int use_magma = 0;

    //int ngpu;

    //ngpu = 0;
    setenv("CUDA_VISIBLE_DEVICES", "3,4", 1);


    IDin = image_ID(IDinput_name);
    IDout = image_ID(IDoutput_name);
    IDRM = image_ID(IDRM_name);


    insize = data.image[IDin].md[0].size[2];
    xsizeout = data.image[IDRM].md[0].size[0];
    ysizeout = data.image[IDRM].md[0].size[1];
    xsizein = data.image[IDin].md[0].size[0];
    ysizein = data.image[IDin].md[0].size[1];

    if(autoMask_MODE == 0)
    {
        IDinmask = image_ID(IDinmask_name);
    }
    else
    {
        create_2Dimage_ID("_RMmask", xsizein, ysizein, &IDinmask);
        for(spl = 0; spl < insize; spl++)
            for(ii = 0; ii < xsizein * ysizein; ii++)
                if(data.image[IDin].array.F[spl * xsizein * ysizein + ii] > 0.5)
                {
                    data.image[IDinmask].array.F[ii] = 1.0;
                }
    }

    // create pokeM
    NBact = 0;
    for(ii = 0; ii < xsizein * ysizein; ii++)
        if(data.image[IDinmask].array.F[ii] > 0.5)
        {
            NBact++;
        }

    printf("NBact = %ld\n", NBact);

    inpixarray = (long *) malloc(sizeof(long) * NBact);
    if(inpixarray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    act = 0;
    for(ii = 0; ii < xsizein * ysizein; ii++)
        if(data.image[IDinmask].array.F[ii] > 0.5)
        {
            inpixarray[act] = ii;
            act++;
        }



    sizearray = (uint32_t *) malloc(sizeof(uint32_t) * 2);
    if(sizearray == NULL) {
        PRINT_ERROR("malloc returns NULL pointer");
        abort();
    }

    sizearray[0] = NBact;
    sizearray[1] = insize; // number of measurements

    printf("NBact = %ld\n", NBact);
    for(act = 0; act < 10; act++)
    {
        printf("act %5ld -> pix %5ld\n", act, inpixarray[act]);
    }


    create_2Dimage_ID("pokeM", NBact, insize, &IDpokeM);

    for(spl = 0; spl < insize; spl++)
        for(act = 0; act < NBact; act++)
        {
            data.image[IDpokeM].array.F[NBact * spl + act] = data.image[IDin].array.F[spl *
                    xsizein * ysizein + inpixarray[act]];
        }
    save_fits("pokeM", "!_test_pokeM.fits");

    // compute pokeM pseudo-inverse
#ifdef HAVE_MAGMA
    CUDACOMP_magma_compute_SVDpseudoInverse("pokeM", "pokeMinv", SVDeps, insize,
                                            "VTmat", 0, 0, 1.e-4, 1.e-7, 0, 64);
#else
    linopt_compute_SVDpseudoInverse("pokeM", "pokeMinv", SVDeps, insize, "VTmat");
#endif

    list_image_ID();
    save_fits("pokeMinv", "!pokeMinv.fits");
    IDpinv = image_ID("pokeMinv");

    // multiply measurements by pokeMinv
    create_3Dimage_ID("_respmat", xsizeout, ysizeout, xsizein * ysizein, &ID_rm);

    for(act = 0; act < NBact; act++)
    {
        for(kk = 0; kk < insize; kk++)
            for(ii = 0; ii < xsizeout * ysizeout; ii++)
            {
                data.image[ID_rm].array.F[inpixarray[act]*xsizeout * ysizeout + ii] +=
                    data.image[IDout].array.F[kk * xsizeout * ysizeout + ii] *
                    data.image[IDpinv].array.F[kk * NBact + act];
            }
    }
    save_fits("_respmat", "!_test_RM.fits");
//exit(0);


    // COMPUTE SOLUTION QUALITY

    IDRM = image_ID("_respmat");


    create_2Dimage_ID("_tmplicli", xsizeout, ysizeout, &IDtmp);
    create_3Dimage_ID("testout", xsizeout, ysizeout, insize, &IDout1);

    printf("IDin  = %ld\n", IDin);
    printf("IDout = %ld\n", IDout);
    printf("IDinmask = %ld\n", IDinmask);

    // on iteration 0, compute initial fit value
    fitval = 0.0;



    for(kk = 0; kk < insize; kk++)
    {
        printf("\r kk = %5ld / %5ld    ", kk, insize);
        fflush(stdout);

        for(ii_out = 0; ii_out < xsizeout; ii_out++)
            for(jj_out = 0; jj_out < ysizeout; jj_out++)
            {
                data.image[IDtmp].array.F[jj_out * xsizeout + ii_out] = 0.0;
            }


        for(ii_in = 0; ii_in < xsizein; ii_in++)
            for(jj_in = 0; jj_in < ysizein; jj_in++)
            {

                //printf("%ld  pix %ld %ld active\n", kk, ii_in, jj_in);
                for(ii_out = 0; ii_out < xsizeout; ii_out++)
                    for(jj_out = 0; jj_out < ysizeout; jj_out++)
                    {
                        data.image[IDtmp].array.F[jj_out * xsizeout + ii_out] +=
                            data.image[IDin].array.F[kk * xsizein * ysizein + jj_in * xsizein + ii_in] *
                            data.image[IDRM].array.F[(jj_in * xsizein + ii_in) * xsizeout * ysizeout +
                                                     jj_out * xsizeout + ii_out];
                    }

            }
        for(ii_out = 0; ii_out < xsizeout; ii_out++)
            for(jj_out = 0; jj_out < ysizeout; jj_out++)
            {
                tmpv1 = data.image[IDtmp].array.F[jj_out * xsizeout + ii_out] -
                        data.image[IDout].array.F[kk * xsizeout * ysizeout + jj_out * xsizeout +
                                                  ii_out];
                fitval += tmpv1 * tmpv1;
                data.image[IDout1].array.F[kk * xsizeout * ysizeout + jj_out * xsizeout +
                                           ii_out] = tmpv1; //data.image[IDtmp].array.F[jj_out*xsizeout+ii_out];
            }
    }
    printf("\n");
    printf("  %5ld    fitval = %.20f\n", kk, sqrt(fitval / xsizeout / ysizeout));

    delete_image_ID("_tmplicli", DELETE_IMAGE_ERRMODE_WARNING);

    free(sizearray);
    free(inpixarray);

    return IDout;
}





























































