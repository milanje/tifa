#if !defined(_TIFA_MAINPAGE_H_)
   /**
    * \def _TIFA_MAINPAGE_H_
    * Standard include guard.
    */
#define _TIFA_MAINPAGE_H_

/**
 * \file    mainpage.h
 * \author  Jerome Milan
 *
 * \mainpage TIFA
 *
 * \section about About the TIFA library 
 * 
 * <strong>TIFA</strong> is an acronym standing for "Tools for Integer
 * FActorization". As its (utterly unoriginal) name implies
 * <strong>TIFA</strong> is an open source library for
 * composite integer factorization. Its goal is to provide portable and 
 * reasonably fast implementations for several algorithms, 
 * with a particular emphasis on the factorization of small to medium-sized 
 * composites, say from 40 bits to about 200 bits. 
 *
 * Although it obviously won't break any record by itself, <strong>TIFA</strong>
 * may be a good companion to more ambitious factorization attempts such as a
 * distributed implementation of the Number Field Sieve, where it can factor the
 * numerous smaller-sized by-products, one per thread: \c tifa_cofactor
 * (cofactor.h) answers whether they split into small enough primes.
 *
 * \section license License
 *
 * Copyright (C) 2011 CNRS - Ecole Polytechnique - INRIA.
 * The files added in 2026 are Copyright (C) 2026 Jerome Milan.
 *
 * This file is part of TIFA.
 *
 * TIFA is free software; you can redistribute it and/or modify it under the 
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version.
 *
 * TIFA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.
 *
 * \section content Content of the TIFA package 
 *
 * Actually, <strong>TIFA</strong> is a little bit more than a library <em>per 
 * se</em>. The <strong>TIFA</strong> package supplies: 
 * 
 * <ul>
 *   <li>a <strong>C99</strong> library providing implementations for the 
 *       following factorization algorithms:
 *     <ul>
 *       <li>CFRAC (Continued FRACtion factorization)</li>
 *       <li>ECM (Elliptic Curve Method)</li>
 *       <li>Fermat (McKee's "fast" variant of Fermat's algorithm)</li>
 *       <li>SIQS (Self-Initializing Quadratic Sieve)</li>
 *       <li>SQUFOF (SQUare FOrm Factorization)</li>
 *       <li>Pollard's rho method</li>
 *       <li>Pollard's p-1 method</li>
 *       <li>Williams's p+1 method</li>
 *       <li>trial division</li>
 *     </ul>
 *     on top of which \c tifa_factor (tifa_factor.h) factors an integer
 *     completely, choosing the methods by size, and \c tifa_cofactor
 *     (cofactor.h) does the cofactorization of sieves. Below 2^64,
 *     primality testing and factorization run in machine words
 *     (factor_u64.h).
 *   </li>
 *   <li>stand-alone factorization programs: <tt>tifa_factor</tt>, which
 *       chooses the methods as the <tt>tifa_factor</tt> function does, and
 *       one program for each algorithm implemented:
 *     <ul>
 *       <li><tt>cfrac</tt></li>
 *       <li><tt>ecm</tt></li>
 *       <li><tt>fermat</tt></li>
 *       <li><tt>siqs</tt></li>
 *       <li><tt>squfof</tt></li>
 *       <li><tt>rho</tt></li>
 *       <li><tt>p_minus_1</tt></li>
 *       <li><tt>p_plus_1</tt></li>
 *       <li><tt>tdiv</tt></li>
 *     </ul>
 *   </li>
 *   <li><tt>tifa_tune</tt>, which measures the algorithms' best parameters on
 *       the current machine.
 *   </li>
 * </ul>
 *
 * \section documentation Documentation
 *
 * Apart from this Doxygen documentation, see the \c README.md file at the top
 * of the source tree, and, in the \c readme directory, \c readme.txt (building,
 * installing, source layout) and \c news.md (the October 2026 rework and its
 * measurements).
 *
 */

#endif
