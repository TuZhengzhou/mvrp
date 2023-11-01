/** @file
 *****************************************************************************

 Declaration of public-parameter selector for vCNN+.

 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef KZG_PARAMS_HPP_
#define KZG_PARAMS_HPP_

#include <libff/algebra/curves/public_params.hpp>

#include <libsnark/relations/constraint_satisfaction_problems/r1cs/r1cs.hpp>

namespace libsnark {

/**
 * Below are various template aliases (used for convenience).
 */

template<typename ppT>
using kzg_constraint_system = r1cs_constraint_system<libff::Fr<ppT> >;

template<typename ppT>
using kzg_primary_input = r1cs_primary_input<libff::Fr<ppT> >;

template<typename ppT>
using kzg_auxiliary_input = r1cs_auxiliary_input<libff::Fr<ppT> >;

} // libsnark

#endif // KZG_PARAMS_HPP_
