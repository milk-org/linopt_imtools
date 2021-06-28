#ifndef LINOPT_IMTOOLS__MASK_TO_PIXTABLE_H
#define LINOPT_IMTOOLS__MASK_TO_PIXTABLE_H

errno_t CLIADDCMD_linopt_imtools__mask_to_pixtable();

long linopt_imtools_mask_to_pixtable(
    const char *IDmask_name,
    const char *IDpixindex_name,
    const char *IDpixmult_name
);

#endif
