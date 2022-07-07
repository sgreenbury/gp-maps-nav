#ifndef _rna_fold_h_included_
#define _rna_fold_h_included_

#ifdef __cplusplus
extern "C" {
#endif

// C declarations
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//#include "utils.h"
#include "include/ViennaRNA/utils.h"
//#include  "fold_vars.h"
#include "include/ViennaRNA/fold.h"
//#include  "part_func.h"
#include "include/ViennaRNA/inverse.h"
//#include  "RNAstruct.h"
//#include  "treedist.h"
//#include  "stringdist.h"
//#include  "profiledist.h"

#ifdef __cplusplus
}
#endif

#include <array>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include "absl/strings/str_format.h"

void Fold(char *&char_seq, char *&phys_seq);
float InverseFold(char *&char_seq, char *&phys_seq);
std::string shape(std::string in, int level);
#endif
