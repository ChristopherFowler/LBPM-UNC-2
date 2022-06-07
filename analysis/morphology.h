// Morphological opening routine
#include "common/Array.h"
#include "common/Domain.h"
#include "analysis/runAnalysis.h"

double MorphOpen(DoubleArray &SignDist, char *id, char* id_original, std::shared_ptr<Domain> Dm, double VoidFraction, char ErodeLabel, char ReplaceLabel, int amin, int amax,
    double deltaR, double Rcrit_new, int count_connected);
