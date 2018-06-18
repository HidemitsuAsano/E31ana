#ifndef KNPIPI_CUT_H
#define KNPIPI_CUT_H 1

#include "GlobalVariables.h"
#include "TMacro.h"

bool isK0(const double im);
bool isIM_Sm(const double im);
bool isIM_Sp(const double im);
bool isSignal(const double pipi, const double npim, const double npip);

void printParam(TMacro *log);

#endif
