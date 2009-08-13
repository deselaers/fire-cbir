#ifndef __kernelfunctionmaker__hpp__
#define __kernelfunctionmaker__hpp__

#include <string>
#include <sstream>

#include "kernelfunction.hpp"
#include "monomial1kernelfunction.hpp"
#include "monomial2kernelfunction.hpp"
#include "monomial3kernelfunction.hpp"
#include "relationalkernelfunction.hpp"

KernelFunction* getNewKernelFunction (const std::string description);

#endif
