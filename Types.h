#pragma once

// LINT is always 64-bit
#define LINT long long int

// useful macros
#define LIMIT(x,xmin,xmax) if (x < xmin) x = xmin; if (x > xmax) x = xmax
#define LIMIT_MIN(x,xmin) if (x < xmin) x = xmin
#define LIMIT_MAX(x,xmax) if (x > xmax) x = xmax
