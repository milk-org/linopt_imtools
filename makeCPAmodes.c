#include <math.h>

// log all debug trace points to file
#define DEBUGLOG

#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_tools/COREMOD_tools.h"

// Local variables pointers
static char *outimname;
static long *sizeout;
static double *CPAmaxval;
static double *deltaCPAval;
static double *radiusval;
static double *radiusfactorlimval;
static long *writefileval;


static CLICMDARGDEF farg[] =
{
    {
        CLIARG_STR, ".out_name", "output image", "out1",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &outimname
    },
    {
        CLIARG_LONG, ".size", "size", "512",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &sizeout
    },
    {
        CLIARG_FLOAT, ".CPAmax", "maximum cycle per aperture", "8.0",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &CPAmaxval
    },
    {
        CLIARG_FLOAT, ".deltaCPA", "CPA interval", "0.8",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &deltaCPAval
    },
    {
        CLIARG_FLOAT, ".radius", "disk radius", "160.0",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &radiusval
    },
    {
        CLIARG_FLOAT, ".radfactlim", "radius factor limit", "1.5",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &radiusfactorlimval
    },
    {
        CLIARG_LONG, ".writefile", "write file flag", "0",
        CLIARG_VISIBLE_DEFAULT,
        (void **) &writefileval
    }
};

static CLICMDDATA CLIcmddata =
{
    "mkFouriermodes",
    "make basis of Fourier Modes",
    CLICMD_FIELDS_DEFAULTS
};


// detailed help
static errno_t help_function()
{
    return RETURN_SUCCESS;
}





errno_t linopt_imtools_makeCPAmodes(
    const char *ID_name,
    long        size,
    float       CPAmax,
    float       deltaCPA,
    float       radius,
    float       radfactlim,
    int         writeMfile,
    long       *outNBmax
)
{
    DEBUG_TRACE_FSTART();
    DEBUG_TRACEPOINT("FARG %s",
                     ID_name);

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
    FUNC_CHECK_RETURN(
        create_2Dimage_ID("cpa_tmpx", size, size, &IDx)
    );

    FUNC_CHECK_RETURN(
        create_2Dimage_ID("cpa_tmpy", size, size, &IDy)
    );

    FUNC_CHECK_RETURN(
        create_2Dimage_ID("cpa_tmpr", size, size, &IDr)
    );

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
    FUNC_CHECK_RETURN(
        create_3Dimage_ID(ID_name, size, size, NBmax - 1, &ID)
    );



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

    FUNC_CHECK_RETURN(
        delete_image_ID("cpamodesfreq", DELETE_IMAGE_ERRMODE_IGNORE)
    );



    DEBUG_TRACEPOINT("Create cpamodesfreq");

    FUNC_CHECK_RETURN(
        create_2Dimage_ID("cpamodesfreq", NBmax - 1, 1, &IDfreq)
    );

    DEBUG_TRACEPOINT("IDfreq %ld", IDfreq);
    DEBUG_TRACEPOINT("IDx %ld", IDx);
    DEBUG_TRACEPOINT("IDy %ld", IDy);
    DEBUG_TRACEPOINT("IDr %ld", IDr);

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
        DEBUG_TRACEPOINT("k = %ld / %ld", k, NBmax);
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
    DEBUG_TRACEPOINT("free memory");

    free(CPAxarray);
    free(CPAyarray);
    free(CPArarray);

    DEBUG_TRACEPOINT("delete tmp files");

    FUNC_CHECK_RETURN(
        delete_image_ID("cpa_tmpx",DELETE_IMAGE_ERRMODE_WARNING)
    );

    FUNC_CHECK_RETURN(
        delete_image_ID("cpa_tmpy", DELETE_IMAGE_ERRMODE_WARNING)
    );

    FUNC_CHECK_RETURN(
        delete_image_ID("cpa_tmpr", DELETE_IMAGE_ERRMODE_WARNING)
    );

    // printf("done \n");
    //fflush(stdout);

    if(outNBmax != NULL)
    {
        *outNBmax = NBmax;
    }

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}




static errno_t compute_function()
{
    DEBUG_TRACE_FSTART();

    INSERT_STD_PROCINFO_COMPUTEFUNC_START

    linopt_imtools_makeCPAmodes(
        outimname,
        *sizeout,
        *CPAmaxval,
        *deltaCPAval,
        *radiusval,
        *radiusfactorlimval,
        *writefileval,
        NULL
    );

    INSERT_STD_PROCINFO_COMPUTEFUNC_END

    DEBUG_TRACE_FEXIT();
    return RETURN_SUCCESS;
}



INSERT_STD_FPSCLIfunctions

// Register function in CLI
errno_t CLIADDCMD_linopt_imtools__makeCPAmodes()
{
    INSERT_STD_CLIREGISTERFUNC
    return RETURN_SUCCESS;
}

