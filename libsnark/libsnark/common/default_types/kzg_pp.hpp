/** @file
 *****************************************************************************

 This file defines default_dario_pp based on the elliptic curve
 choice selected in ec_pp.hpp.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef KZG_PP_HPP_
#define KZG_PP_HPP_

#include <libff/common/default_types/ec_pp.hpp>

namespace libsnark {
typedef libff::default_ec_pp default_kzg_pp;
} // libsnark

#endif // KZG_PP_HPP_
