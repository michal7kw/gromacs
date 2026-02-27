# Changes to the API

## Legacy aggregating headers have been removed.

Previously, some of the legacy API headers existed
only\_ to aggregate `#include` lines for other
installed headers. No guidance was provided regarding which header to
include for a given feature. These redundant headers have been removed.
Client software relying on `#include "gromacs/module.h"` will need to be
updated with more specific `#include "gromacs/module/feature.h"`
directives.

`4487`
