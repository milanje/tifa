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
 * FActorisation". As its (utterly unoriginal) name implies 
 * <strong>TIFA</strong> is a open source library for
 * composite integer factorization. Its goal is to provide portable and 
 * reasonably fast implementations for several algorithms, 
 * with a particular emphasis on the factorization of small to medium-sized 
 * composites, say from 40 bits to about 200 bits. 
 *
 * Although it obviously won't break any record by itself, <strong>TIFA</strong> 
 * may be a good companion to more ambitious factorization attempts such as a 
 * distributed implementation of the Number Field Sieve, where it could be used 
 * to factor the numerous smaller-sized by-products. 
 *
 * \section license License
 *
 * Copyright (C) 2011 CNRS - Ecole Polytechnique - INRIA.
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
 *     </ul>
 *   </li>
 *   <li>a set of stand-alone factorization programs for each algorithm
 *       implemented:
 *     <ul>
 *       <li><tt>cfrac_program</tt></li>
 *       <li><tt>ecm_program</tt></li>
 *       <li><tt>fermat_program</tt></li>
 *       <li><tt>siqs_program</tt></li> 
 *       <li><tt>squfof_program</tt></li>
 *     </ul> 
 *   </li>
 *   <li>a set of <strong>Perl 5</strong> scripts wrappers and launchers;</li>
 *   <li>a basic benchmarking framework written in <strong>Perl 5</strong> used 
 *       to assess the performance of <strong>TIFA</strong>'s implementations.
 *   </li> 
 * </ul> 
 *
 * \section documentation Documentation
 *
 * A complete user's guide is in preparation.
 *
 * In the interim, the best source of documentation (apart from this Doxygen
 * documentation generated during the build process) is the included (infamous)
 * \c readme.txt file in the \c readme directory.
 *
 * Also worth a look is the (unfortunately not empty) \c issues.txt file.
 *
 */

#endif
