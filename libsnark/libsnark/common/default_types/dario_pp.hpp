/** @file
 *****************************************************************************

 This file defines default_dario_pp based on the elliptic curve
 choice selected in ec_pp.hpp.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef DARIO_PP_HPP_
#define DARIO_PP_HPP_

#include <libff/common/default_types/ec_pp.hpp>

namespace libsnark {
typedef libff::default_ec_pp default_dario_pp;
} // libsnark

#endif // DARIO_PP_HPP_
