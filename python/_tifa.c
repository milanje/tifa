//
// Copyright (C) 2026 Jerome Milan.
//
// This file is part of TIFA.
//
// TIFA is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// TIFA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//

//
// tifa._tifa: tifa_factor_with_hint, tifa_prove_prime, tifa_batch_gcd,
// tifa_cofactor, tifa_tuning_set_dir, tifa_set_time_limit and
// tifa_set_interrupt_check for Python. The argument
// checks and the shape of the results are left to tifa/__init__.py.
//
#if !defined(Py_LIMITED_API)
#define Py_LIMITED_API 0x030A0000
#endif
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <limits.h>
#include <string.h>
#include <time.h>

#include <gmp.h>

#include "tifa_config.h"
#include "array.h"
#include "batch_gcd.h"
#include "cofactor.h"
#include "exit_codes.h"
#include "factoring_machine.h"
#include "prime_cert.h"
#include "tifa_factor.h"
#include "time_limit.h"
#include "tuning.h"

//
// The signal handlers (KeyboardInterrupt on Ctrl-C) run in the main thread
// when PyErr_CheckSignals is called, which needs the GIL: the factorizations
// ask at most this often, their thread having released it.
//
#define SIGNAL_POLL_SECONDS 0.05

static double monotonic_now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double) t.tv_sec + 1e-9 * (double) t.tv_nsec;
}

static bool check_signals(void* const arg) {
    double* const next = arg;
    const double now = monotonic_now();
    if (now < *next) {
        return false;
    }
    *next = now + SIGNAL_POLL_SECONDS;
    const PyGILState_STATE gil = PyGILState_Ensure();
    const int raised = PyErr_CheckSignals();
    PyGILState_Release(gil);
    return raised != 0;
}

//
// Integers cross the boundary as hexadecimal strings, a few hundred
// nanoseconds, unless they fit in a long: a cofactorization can take a
// microsecond.
//
static int mpz_from_pylong(mpz_t z, PyObject* const n) {
    int overflow;
    const long long value = PyLong_AsLongLongAndOverflow(n, &overflow);
    if ((value == -1) && PyErr_Occurred()) {
        return -1;
    }
    if ((overflow == 0) && (value >= LONG_MIN) && (value <= LONG_MAX)) {
        mpz_set_si(z, (long) value);
        return 0;
    }
    PyObject* const hex = PyNumber_ToBase(n, 16);
    if (hex == NULL) {
        return -1;
    }
    const char* const digits = PyUnicode_AsUTF8AndSize(hex, NULL);
    const int ok = (digits != NULL) && (mpz_set_str(z, digits, 0) == 0);
    Py_DECREF(hex);
    if (!ok && !PyErr_Occurred()) {
        PyErr_SetString(PyExc_ValueError, "not an integer");
    }
    return ok ? 0 : -1;
}

static PyObject* pylong_from_mpz(const mpz_t z) {
    if (mpz_fits_slong_p(z)) {
        return PyLong_FromLong(mpz_get_si(z));
    }
    char* const digits = PyMem_Malloc(mpz_sizeinbase(z, 16) + 2);
    if (digits == NULL) {
        return PyErr_NoMemory();
    }
    mpz_get_str(digits, 16, z);
    PyObject* const n = PyLong_FromString(digits, NULL, 16);
    PyMem_Free(digits);
    return n;
}

//
// The ints of a list, all positive (what = "moduli": "moduli must be ...").
//
static mpz_array_t* positive_ints(PyObject* const list,
                                  const char* const what) {
    const Py_ssize_t count = PyList_Size(list);
    if (count > (Py_ssize_t) UINT32_MAX - ELONGATION) {
        PyErr_Format(PyExc_OverflowError, "too many %s", what);
        return NULL;
    }
    mpz_array_t* const ints = alloc_mpz_array((uint32_t) count);
    for (Py_ssize_t i = 0; i < count; i++) {
        PyObject* const n = PyList_GetItem(list, i);
        if (   !PyLong_Check(n)
            || (mpz_from_pylong(ints->data[i], n) != 0)
            || (mpz_sgn(ints->data[i]) <= 0)) {
            if (!PyErr_Occurred()) {
                PyErr_Format(PyExc_ValueError, "%s must be positive ints",
                             what);
            }
            free_mpz_array(ints);
            return NULL;
        }
        ints->length++;
    }
    return ints;
}

static PyObject* pairs_from_arrays(const mpz_array_t* const factors,
                                   const uint32_array_t* const multis) {
    PyObject* const pairs = PyList_New(0);
    if (pairs == NULL) {
        return NULL;
    }
    for (uint32_t i = 0; i < factors->length; i++) {
        PyObject* const p = pylong_from_mpz(factors->data[i]);
        PyObject* const pair = (p == NULL)
                             ? NULL
                             : Py_BuildValue("(NI)", p, multis->data[i]);
        if ((pair == NULL) || (PyList_Append(pairs, pair) != 0)) {
            Py_XDECREF(pair);
            Py_DECREF(pairs);
            return NULL;
        }
        Py_DECREF(pair);
    }
    return pairs;
}

PyDoc_STRVAR(factor_doc,
"_factor(n, hint, timeout, interruptible, /)\n--\n\n"
"Factor the int n > 1 with tifa_factor_with_hint (hint 0, 1 or 2), within\n"
"timeout seconds if it is positive and finite, and stopping for the\n"
"exception of a signal handler if interruptible (in the main thread). Return\n"
"(complete, [(factor, exponent), ...], timed_out), complete being False when\n"
"some factors may not be prime.");

static PyObject* factor(PyObject* const self, PyObject* const args) {
    (void) self;
    PyObject* nobj;
    int hint;
    double timeout;
    int interruptible;
    if (!PyArg_ParseTuple(args, "O!idp:_factor", &PyLong_Type, &nobj, &hint,
                          &timeout, &interruptible)) {
        return NULL;
    }
    mpz_t n;
    mpz_init(n);
    if (mpz_from_pylong(n, nobj) != 0) {
        mpz_clear(n);
        return NULL;
    }
    if ((mpz_cmp_ui(n, 1) <= 0) || (hint < 0) || (hint > 2)) {
        mpz_clear(n);
        PyErr_SetString(PyExc_ValueError, "n must be above 1, hint in [0, 2]");
        return NULL;
    }
    mpz_array_t*    const factors = alloc_mpz_array(8);
    uint32_array_t* const multis  = alloc_uint32_array(8);
    ecode_t ecode;
    int prime = 0;
    bool timed_out;
    bool interrupted;
    double next_poll = monotonic_now() + SIGNAL_POLL_SECONDS;

    Py_BEGIN_ALLOW_THREADS
    tifa_set_time_limit(timeout);
    if (interruptible) {
        tifa_set_interrupt_check(check_signals, &next_poll);
    }
    ecode = tifa_factor_with_hint(factors, multis, n,
                                  FIND_COMPLETE_FACTORIZATION,
                                  (tifa_factor_hint_t) hint);
    timed_out   = tifa_time_limit_reached();
    interrupted = tifa_interrupted();
    tifa_set_interrupt_check(NULL, NULL);
    tifa_set_time_limit(0.0);
    if ((ecode == NO_FACTOR_FOUND) && !timed_out && !interrupted) {
        prime = (mpz_probab_prime_p(n, 30) != 0);
    }
    Py_END_ALLOW_THREADS

    if (interrupted) {
        free_mpz_array(factors);
        free_uint32_array(multis);
        mpz_clear(n);
        return NULL;
    }

    PyObject* const late = timed_out ? Py_True : Py_False;
    PyObject* result = NULL;
    if (prime) {
        PyObject* const p = pylong_from_mpz(n);
        if (p != NULL) {
            result = Py_BuildValue("(O[(NI)]O)", Py_True, p, 1U, late);
        }
    } else if (   (ecode == COMPLETE_FACTORIZATION_FOUND)
               || (ecode == PARTIAL_FACTORIZATION_FOUND)) {
        PyObject* const pairs = pairs_from_arrays(factors, multis);
        if (pairs != NULL) {
            result = Py_BuildValue(
                "(ONO)",
                (ecode == COMPLETE_FACTORIZATION_FOUND) ? Py_True : Py_False,
                pairs,
                late
            );
        }
    } else {
        PyObject* const p = pylong_from_mpz(n);
        if (p != NULL) {
            result = Py_BuildValue("(O[(NI)]O)", Py_False, p, 1U, late);
        }
    }
    free_mpz_array(factors);
    free_uint32_array(multis);
    mpz_clear(n);
    return result;
}

static PyObject* witness_pairs(const mpz_array_t* const primes,
                               const uint32_array_t* const witnesses) {
    PyObject* const pairs = PyList_New(0);
    if (pairs == NULL) {
        return NULL;
    }
    for (uint32_t j = 0; j < primes->length; j++) {
        PyObject* const q = pylong_from_mpz(primes->data[j]);
        PyObject* const pair = (q == NULL)
                             ? NULL
                             : Py_BuildValue("(NI)", q, witnesses->data[j]);
        if ((pair == NULL) || (PyList_Append(pairs, pair) != 0)) {
            Py_XDECREF(pair);
            Py_DECREF(pairs);
            return NULL;
        }
        Py_DECREF(pair);
    }
    return pairs;
}

static PyObject* step_from_cert(const prime_cert_step_t* const step) {
    PyObject* const p = pylong_from_mpz(step->p);
    PyObject* const qs = witness_pairs(step->pm1_primes, step->pm1_witnesses);
    PyObject* const rs = (step->pp1_d == 0)
                       ? NULL
                       : witness_pairs(step->pp1_primes, step->pp1_witnesses);
    if ((p == NULL) || (qs == NULL) || ((step->pp1_d != 0) && (rs == NULL))) {
        Py_XDECREF(p);
        Py_XDECREF(qs);
        Py_XDECREF(rs);
        return NULL;
    }
    return (step->pp1_d == 0)
         ? Py_BuildValue("(NN)", p, qs)
         : Py_BuildValue("(NNiN)", p, qs, (int) step->pp1_d, rs);
}

PyDoc_STRVAR(prime_certificate_doc,
"_prime_certificate(n, timeout, interruptible, /)\n--\n\n"
"Prove the int n prime with tifa_prove_prime, within timeout seconds if it is\n"
"positive and finite, and stopping for the exception of a signal handler if\n"
"interruptible. Return (status, steps, timed_out), status 0 (proven), 1 (not\n"
"prime) or 2 (not proven), steps [(p, [(q, a), ...]), ...] or\n"
"(p, [(q, a), ...], D, [(r, P), ...]) when proven.");

static PyObject* prime_certificate(PyObject* const self, PyObject* const args) {
    (void) self;
    PyObject* nobj;
    double timeout;
    int interruptible;
    if (!PyArg_ParseTuple(args, "O!dp:_prime_certificate", &PyLong_Type,
                          &nobj, &timeout, &interruptible)) {
        return NULL;
    }
    mpz_t n;
    mpz_init(n);
    if (mpz_from_pylong(n, nobj) != 0) {
        mpz_clear(n);
        return NULL;
    }
    prime_cert_t* const cert = alloc_prime_cert();
    prime_cert_status_t status;
    bool timed_out;
    bool interrupted;
    double next_poll = monotonic_now() + SIGNAL_POLL_SECONDS;

    Py_BEGIN_ALLOW_THREADS
    tifa_set_time_limit(timeout);
    if (interruptible) {
        tifa_set_interrupt_check(check_signals, &next_poll);
    }
    status      = tifa_prove_prime(cert, n);
    timed_out   = tifa_time_limit_reached();
    interrupted = tifa_interrupted();
    tifa_set_interrupt_check(NULL, NULL);
    tifa_set_time_limit(0.0);
    Py_END_ALLOW_THREADS

    if (interrupted) {
        free_prime_cert(cert);
        mpz_clear(n);
        return NULL;
    }

    PyObject* steps = PyList_New(0);
    for (uint32_t i = 0; (steps != NULL) && (i < cert->nsteps); i++) {
        PyObject* const step = step_from_cert(&cert->steps[i]);
        if ((step == NULL) || (PyList_Append(steps, step) != 0)) {
            Py_CLEAR(steps);
        }
        Py_XDECREF(step);
    }
    free_prime_cert(cert);
    mpz_clear(n);
    return (steps == NULL) ? NULL
                           : Py_BuildValue("(iNO)", (int) status, steps,
                                           timed_out ? Py_True : Py_False);
}

PyDoc_STRVAR(batch_gcd_doc,
"_batch_gcd(moduli, timeout, interruptible, /)\n--\n\n"
"tifa_batch_gcd on the list of ints moduli, all positive, within timeout\n"
"seconds if it is positive and finite, and stopping for the exception of a\n"
"signal handler if interruptible. Return ([gcd, ...], timed_out), the list\n"
"empty if timed out.");

static PyObject* batch_gcd(PyObject* const self, PyObject* const args) {
    (void) self;
    PyObject* list;
    double timeout;
    int interruptible;
    if (!PyArg_ParseTuple(args, "O!dp:_batch_gcd", &PyList_Type, &list,
                          &timeout, &interruptible)) {
        return NULL;
    }
    mpz_array_t* const moduli = positive_ints(list, "moduli");
    if (moduli == NULL) {
        return NULL;
    }
    mpz_array_t* const gcds = alloc_mpz_array(moduli->length);
    ecode_t ecode;
    bool timed_out;
    bool interrupted;
    double next_poll = monotonic_now() + SIGNAL_POLL_SECONDS;

    Py_BEGIN_ALLOW_THREADS
    tifa_set_time_limit(timeout);
    if (interruptible) {
        tifa_set_interrupt_check(check_signals, &next_poll);
    }
    ecode       = tifa_batch_gcd(gcds, moduli);
    timed_out   = tifa_time_limit_reached();
    interrupted = tifa_interrupted();
    tifa_set_interrupt_check(NULL, NULL);
    tifa_set_time_limit(0.0);
    Py_END_ALLOW_THREADS

    free_mpz_array(moduli);
    if (interrupted) {
        free_mpz_array(gcds);
        return NULL;
    }
    PyObject* result = PyList_New(0);
    for (uint32_t i = 0; (result != NULL) && (i < gcds->length); i++) {
        PyObject* const g = pylong_from_mpz(gcds->data[i]);
        if ((g == NULL) || (PyList_Append(result, g) != 0)) {
            Py_CLEAR(result);
        }
        Py_XDECREF(g);
    }
    free_mpz_array(gcds);
    if ((result != NULL) && (ecode == FAILURE)) {
        Py_DECREF(result);
        PyErr_SetString(PyExc_ValueError, "moduli must be positive ints");
        return NULL;
    }
    return (result == NULL) ? NULL
                            : Py_BuildValue("(NO)", result,
                                            timed_out ? Py_True : Py_False);
}

#define STRATEGY_NAME "tifa._tifa.cofactor_strategy"

static void free_strategy(PyObject* const capsule) {
    free_cofactor_strategy(PyCapsule_GetPointer(capsule, STRATEGY_NAME));
}

PyDoc_STRVAR(cofactor_strategy_doc,
"_cofactor_strategy(lpb, fb_bound, max_large_primes, max_siqs_bits, /)\n--\n\n"
"alloc_cofactor_strategy for these bounds (unchecked but lpb), in a capsule.");

static PyObject* cofactor_strategy(PyObject* const self, PyObject* const args) {
    (void) self;
    cofactor_params_t params;
    unsigned long long fb_bound;
    if (!PyArg_ParseTuple(args, "IKII:_cofactor_strategy", &params.lpb,
                          &fb_bound, &params.max_large_primes,
                          &params.max_siqs_bits)) {
        return NULL;
    }
    params.fb_bound = (uint64_t) fb_bound;
    cofactor_strategy_t* strategy;
    Py_BEGIN_ALLOW_THREADS
    strategy = alloc_cofactor_strategy(&params);
    Py_END_ALLOW_THREADS
    if (strategy == NULL) {
        PyErr_SetString(PyExc_ValueError, "lpb must be in [1, 63]");
        return NULL;
    }
    PyObject* const capsule = PyCapsule_New(strategy, STRATEGY_NAME,
                                            free_strategy);
    if (capsule == NULL) {
        free_cofactor_strategy(strategy);
    }
    return capsule;
}

//
// The primes of all the smooth numbers of a call, appended in turn: their
// number per smooth number is kept apart.
//
typedef struct {
    mpz_array_t*    primes;
    uint32_array_t* exponents;
} found_t;

static void keep_found(found_t* const found, const mpz_array_t* const primes,
                       const uint32_array_t* const exponents) {
    const uint32_t needed = found->primes->length + primes->length;
    if (needed > found->primes->alloced) {
        const uint32_t alloced = (needed > UINT32_MAX / 2) ? needed
                                                           : 2 * needed;
        resize_mpz_array(found->primes, alloced);
        resize_uint32_array(found->exponents, alloced);
    }
    for (uint32_t k = 0; k < primes->length; k++) {
        mpz_set(found->primes->data[found->primes->length++],
                primes->data[k]);
        found->exponents->data[found->exponents->length++] =
            exponents->data[k];
    }
}

//
// {prime: exponent} of the count primes from first, in increasing order.
//
static PyObject* factor_dict(const found_t* const found, const uint32_t first,
                             const uint32_t count) {
    uint32_t order[64];
    uint32_t* const index = (count <= 64)
                          ? order : PyMem_Malloc(count * sizeof(uint32_t));
    if (index == NULL) {
        return PyErr_NoMemory();
    }
    for (uint32_t k = 0; k < count; k++) {
        uint32_t j = k;
        while (   (j > 0)
               && (mpz_cmp(found->primes->data[first + index[j - 1]],
                           found->primes->data[first + k]) > 0)) {
            index[j] = index[j - 1];
            j--;
        }
        index[j] = k;
    }
    PyObject* dict = PyDict_New();
    for (uint32_t k = 0; (dict != NULL) && (k < count); k++) {
        PyObject* const p = pylong_from_mpz(
                                found->primes->data[first + index[k]]
                            );
        PyObject* const e = (p == NULL)
                          ? NULL
                          : PyLong_FromUnsignedLong(
                                found->exponents->data[first + index[k]]
                            );
        if ((e == NULL) || (PyDict_SetItem(dict, p, e) != 0)) {
            Py_CLEAR(dict);
        }
        Py_XDECREF(p);
        Py_XDECREF(e);
    }
    if (index != order) {
        PyMem_Free(index);
    }
    return dict;
}

PyDoc_STRVAR(cofactor_doc,
"_cofactor(strategy, ns, interruptible, /)\n--\n\n"
"tifa_cofactor on each int of the list ns, all positive, with the strategy\n"
"of _cofactor_strategy, stopping for the exception of a signal handler if\n"
"interruptible. Return a list of {prime: exponent} (primes in increasing\n"
"order) for the smooth ints, None for the others.");

static PyObject* cofactor(PyObject* const self, PyObject* const args) {
    (void) self;
    PyObject* capsule;
    PyObject* list;
    int interruptible;
    if (!PyArg_ParseTuple(args, "OO!p:_cofactor", &capsule, &PyList_Type,
                          &list, &interruptible)) {
        return NULL;
    }
    const cofactor_strategy_t* const strategy =
        PyCapsule_GetPointer(capsule, STRATEGY_NAME);
    if (strategy == NULL) {
        return NULL;
    }
    mpz_array_t* const ns = positive_ints(list, "ns");
    if (ns == NULL) {
        return NULL;
    }
    const uint32_t count = ns->length;
    //
    // Number of primes of each n, or UINT32_MAX if it is not smooth.
    //
    uint32_t* const nprimes = PyMem_Malloc((count + 1) * sizeof(uint32_t));
    if (nprimes == NULL) {
        free_mpz_array(ns);
        return PyErr_NoMemory();
    }
    found_t found;
    found.primes    = alloc_mpz_array(ELONGATION);
    found.exponents = alloc_uint32_array(ELONGATION);
    mpz_array_t*    const primes    = alloc_mpz_array(ELONGATION);
    uint32_array_t* const exponents = alloc_uint32_array(ELONGATION);
    bool interrupted = false;
    double next_poll = monotonic_now() + SIGNAL_POLL_SECONDS;

    Py_BEGIN_ALLOW_THREADS
    if (interruptible) {
        tifa_set_interrupt_check(check_signals, &next_poll);
    }
    for (uint32_t i = 0; i < count; i++) {
        nprimes[i] = UINT32_MAX;
        if (tifa_cofactor(primes, exponents, ns->data[i], strategy)
                == COFACTOR_SMOOTH) {
            keep_found(&found, primes, exponents);
            nprimes[i] = primes->length;
        }
        if (   tifa_interrupted()
            || (interruptible && check_signals(&next_poll))) {
            interrupted = true;
            break;
        }
    }
    tifa_set_interrupt_check(NULL, NULL);
    Py_END_ALLOW_THREADS

    free_mpz_array(primes);
    free_uint32_array(exponents);
    free_mpz_array(ns);

    PyObject* result = interrupted ? NULL : PyList_New(count);
    uint32_t first = 0;
    for (uint32_t i = 0; (result != NULL) && (i < count); i++) {
        PyObject* item;
        if (nprimes[i] == UINT32_MAX) {
            item = Py_None;
            Py_INCREF(item);
        } else {
            item = factor_dict(&found, first, nprimes[i]);
            first += nprimes[i];
        }
        if ((item == NULL) || (PyList_SetItem(result, i, item) != 0)) {
            Py_CLEAR(result);
        }
    }
    PyMem_Free(nprimes);
    free_mpz_array(found.primes);
    free_uint32_array(found.exponents);
    return result;
}

PyDoc_STRVAR(set_tuning_dir_doc,
"_set_tuning_dir(path, /)\n--\n\n"
"tifa_tuning_set_dir(path), path being bytes without NUL, or None.");

static PyObject* set_tuning_dir(PyObject* const self, PyObject* const path) {
    (void) self;
    if (path == Py_None) {
        tifa_tuning_set_dir(NULL);
        Py_RETURN_NONE;
    }
    char* dir;
    Py_ssize_t length;
    if (PyBytes_AsStringAndSize(path, &dir, &length) != 0) {
        return NULL;
    }
    if (strlen(dir) != (size_t) length) {
        PyErr_SetString(PyExc_ValueError, "embedded null byte");
        return NULL;
    }
    tifa_tuning_set_dir(dir);
    Py_RETURN_NONE;
}

static PyMethodDef methods[] = {
    {"_factor", factor, METH_VARARGS, factor_doc},
    {"_prime_certificate", prime_certificate, METH_VARARGS,
     prime_certificate_doc},
    {"_set_tuning_dir", set_tuning_dir, METH_O, set_tuning_dir_doc},
    {"_batch_gcd", batch_gcd, METH_VARARGS, batch_gcd_doc},
    {"_cofactor_strategy", cofactor_strategy, METH_VARARGS,
     cofactor_strategy_doc},
    {"_cofactor", cofactor, METH_VARARGS, cofactor_doc},
    {NULL, NULL, 0, NULL}
};

static int exec_module(PyObject* const module) {
    return PyModule_AddStringConstant(module, "LIBRARY_VERSION", TIFA_VERSION);
}

static PyModuleDef_Slot slots[] = {
    {Py_mod_exec, (void*) exec_module},
    {0, NULL}
};

static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "_tifa",
    "The C part of the tifa package.",
    0,
    methods,
    slots,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit__tifa(void) {
    return PyModuleDef_Init(&module_def);
}
