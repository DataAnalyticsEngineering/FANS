/*
 * C interface of a FANS material plugin: a library FANS dlopen()s at run time
 * (libfans_<name>.so for "matmodel": "<name>"), so the plugin's dependencies
 * -- libtorch, for NEML2 -- never reach the FANS build.
 *
 * Strain and stress are [n_points][6] in FANS Mandel order,
 * [11, 22, 33, sqrt(2)*12, sqrt(2)*13, sqrt(2)*23].
 * A model may also take a crystal orientation, [n_points][3][3] row-major
 * crystal-to-sample rotation matrices, which FANS reads per grain from the
 * "rotation_matrices" dataset next to the microstructure.
 * History is opaque: n_state doubles per point that FANS stores, starts at
 * zero, and hands back as state_old once a step has converged.
 * Functions returning int give 0 on success, else a message in `err`.
 */

#ifndef FANS_PLUGIN_H
#define FANS_PLUGIN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define FANS_PLUGIN_ERRLEN 1024

typedef struct FANSPluginModel FANSPluginModel;

/* `spec` is the "artifact" material property. NULL on failure. */
FANSPluginModel *fans_plugin_load(const char *spec, int *n_state, int *wants_orientation,
                                  char *err, size_t errlen);

/* orientation is NULL unless wanted; state_old/state_new are NULL when n_state == 0. */
int fans_plugin_evaluate(FANSPluginModel *model, size_t n_points, const double *strain,
                         const double *orientation, const double *state_old, double *stress,
                         double *state_new, char *err, size_t errlen);

void fans_plugin_free(FANSPluginModel *model);

#ifdef __cplusplus
}
#endif
#endif /* FANS_PLUGIN_H */
