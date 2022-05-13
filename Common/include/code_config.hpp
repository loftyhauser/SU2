/*!
 * \file code_config.hpp
 * \brief Header file for collecting common macros, definitions and type configurations.
 * \author T. Albring, P. Gomes, J. Blühdorn
 * \version 7.3.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <type_traits>

#if defined(_MSC_VER)
#define FORCEINLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define FORCEINLINE inline __attribute__((always_inline))
#else
#define FORCEINLINE inline
#endif

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define NEVERINLINE inline __attribute__((noinline))
#else
#define NEVERINLINE inline
#endif

#if defined(__INTEL_COMPILER)
/*--- Disable warnings related to inline attributes. ---*/
#pragma warning disable 2196
#pragma warning disable 3415
/*--- Disable warnings related to overloaded virtual. ---*/
#pragma warning disable 654
#pragma warning disable 1125
#endif

/*--- Convenience SFINAE typedef to conditionally
 * enable/disable function template overloads. ---*/
template<bool condition>
using su2enable_if = typename std::enable_if<condition,bool>::type;

/*--- Compile-time type selection. ---*/
template<bool B, class T, class F> struct su2conditional { using type = T; };
template<class T, class F> struct su2conditional<false, T, F> { using type = F; };

template<bool B, class T, class F>
using su2conditional_t = typename su2conditional<B,T,F>::type;

/*! \brief Static cast "In" to "Out", in debug builds a dynamic cast is used. */
template<class Out, class In>
FORCEINLINE Out su2staticcast_p(In ptr) {
  static_assert(std::is_pointer<In>::value, "This expects a pointer");
#ifndef NDEBUG
  return static_cast<Out>(ptr);
#else
  return dynamic_cast<Out>(ptr);
#endif
}

/*--- Detect compilation with OpenMP. ---*/
#if defined(_OPENMP)
#define HAVE_OMP
#endif

/*--- Depending on the datatype defined during the configuration,
 * include the correct definition, and create the main typedef. ---*/

using su2double = double;

/*--- This type can be used for (rare) compatibility cases or for
 * computations that are intended to be (always) passive. ---*/
using passivedouble = double;

/*--- Define a type for potentially lower precision operations. ---*/
using su2mixedfloat = passivedouble;
