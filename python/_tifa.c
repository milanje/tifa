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
// tifa._tifa: tifa_factor_with_hint, tifa_prove_prime and tifa_tuning_set_dir
// for Python. The argument checks and the shape of the results are left to
// tifa/__init__.py.
//
#if !defined(Py_LIMITED_API)
#define Py_LIMITED_API 0x030A0000
#endif
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <string.h>

#include <gmp.h>

#include "tifa_config.h"
#include "array.h"
#include "exit_codes.h"
#include "factoring_machine.h"
#include "prime_cert.h"
#include "tifa_factor.h"
#include "tuning.h"

//
// Integers cross the boundary as hexadecimal strings: a few hundred
// nanoseconds, against microseconds for the smallest factorization.
//
static int mpz_from_pylong(mpz_t z, PyObject* const n) {
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
    char* const digits = PyMem_Malloc(mpz_sizeinbase(z, 16) + 2);
    if (digits == NULL) {
        return PyErr_NoMemory();
    }
    mpz_get_str(digits, 16, z);
    PyObject* const n = PyLong_FromString(digits, NULL, 16);
    PyMem_Free(digits);
    return n;
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
"_factor(n, hint, /)\n--\n\n"
"Factor the int n > 1 with tifa_factor_with_hint (hint 0, 1 or 2).\n"
"Return (complete, [(factor, exponent), ...]), complete being False when\n"
"some factors may not be prime.");

static PyObject* factor(PyObject* const self, PyObject* const args) {
    (void) self;
    PyObject* nobj;
    int hint;
    if (!PyArg_ParseTuple(args, "O!i:_factor", &PyLong_Type, &nobj, &hint)) {
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

    Py_BEGIN_ALLOW_THREADS
    ecode = tifa_factor_with_hint(factors, multis, n,
                                  FIND_COMPLETE_FACTORIZATION,
                                  (tifa_factor_hint_t) hint);
    if (ecode == NO_FACTOR_FOUND) {
        prime = (mpz_probab_prime_p(n, 30) != 0);
    }
    Py_END_ALLOW_THREADS

    PyObject* result = NULL;
    if (prime) {
        PyObject* const p = pylong_from_mpz(n);
        if (p != NULL) {
            result = Py_BuildValue("(O[(NI)])", Py_True, p, 1U);
        }
    } else if (   (ecode == COMPLETE_FACTORIZATION_FOUND)
               || (ecode == PARTIAL_FACTORIZATION_FOUND)) {
        PyObject* const pairs = pairs_from_arrays(factors, multis);
        if (pairs != NULL) {
            result = Py_BuildValue(
                "(ON)",
                (ecode == COMPLETE_FACTORIZATION_FOUND) ? Py_True : Py_False,
                pairs
            );
        }
    } else {
        PyObject* const p = pylong_from_mpz(n);
        if (p != NULL) {
            result = Py_BuildValue("(O[(NI)])", Py_False, p, 1U);
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
"_prime_certificate(n, /)\n--\n\n"
"Prove the int n prime with tifa_prove_prime. Return (status, steps), status\n"
"0 (proven), 1 (not prime) or 2 (not proven), steps [(p, [(q, a), ...]), ...]\n"
"or (p, [(q, a), ...], D, [(r, P), ...]) when proven.");

static PyObject* prime_certificate(PyObject* const self, PyObject* const args) {
    (void) self;
    PyObject* nobj;
    if (!PyArg_ParseTuple(args, "O!:_prime_certificate", &PyLong_Type,
                          &nobj)) {
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

    Py_BEGIN_ALLOW_THREADS
    status = tifa_prove_prime(cert, n);
    Py_END_ALLOW_THREADS

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
    return (steps == NULL) ? NULL : Py_BuildValue("(iN)", (int) status, steps);
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
