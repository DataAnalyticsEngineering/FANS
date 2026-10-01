/*
 * C interface of a FANS material plugin: a library FANS dlopen()s at run time
 * (libfans_<name>.so for "matmodel": "<name>"), so the plugin's dependencies
 * -- libtorch, for NEML2 -- never reach the FANS build.
 *
 * A plugin material maps a gradient to a flux: the strain to the stress
 * (small strain, Mandel order [11, 22, 33, sqrt(2)*12, sqrt(2)*13, sqrt(2)*23]),
 * the deformation gradient F to the first Piola stress P (large strain,
 * row-major 3x3), or the temperature gradient to the heat flux (thermal).
 * All arrays are in host memory, whatever device the plugin computes on, and
 * hold one row per point:
 *   gradient, flux  [n_points][n_str]
 *   fields[f]       [n_points][size of field f]: further model inputs, e.g. an
 *                   orientation, which FANS reads from the dataset of their
 *                   name next to the microstructure, per voxel [Z][Y][X][...]
 *                   or per phase [n_phase][...]; the "fields" material property
 *                   may name another dataset, {"orientation": "rotation_matrices"}
 *   history_old/new [n_points][sum of the history sizes]: internal variables
 *                   that FANS stores, starts at zero, and hands back as
 *                   history_old once a step has converged; history_new arrives
 *                   holding the latest trial values, e.g. as an initial guess
 * Symmetric tensors among the fields and history are in FANS Mandel order too.
 * Functions returning int give 0 on success, else a message in `err`.
 */

#ifndef FANS_PLUGIN_H
#define FANS_PLUGIN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define FANS_PLUGIN_MSGLEN 1024

typedef struct FANSPluginModel FANSPluginModel;

/* `config` is the material's "material_properties" as JSON. On success `msg`
   gets, as JSON, the gradient's size, the points per fans_plugin_evaluate,
   and the history variables and fields in storage order as [name, size], e.g.
   {"gradient": 6, "batch_size": 65536, "history": [["state/internal/Ep", 6]],
    "fields": [["orientation", 9]]}.
   On failure (NULL) it gets the error. */
FANSPluginModel *fans_plugin_load(const char *config, char *msg, size_t msglen);

/* The step runs from time t_old to t. */
int fans_plugin_evaluate(FANSPluginModel *model, size_t n_points, double t_old, double t, const double *gradient,
                         const double *const *fields, const double *history_old, double *flux,
                         double *history_new, char *err, size_t errlen);

void fans_plugin_free(FANSPluginModel *model);

#ifdef __cplusplus
}
#endif
#endif /* FANS_PLUGIN_H */
