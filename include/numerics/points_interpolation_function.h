// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef LIBMESH_POINTS_INTERPOLATION_FUNCTION_H
#define LIBMESH_POINTS_INTERPOLATION_FUNCTION_H

#include "libmesh/libmesh_config.h"
#include "libmesh/function_base.h"

// Local includes
#include "libmesh/numeric_vector.h"
#include "libmesh/threads.h"
#include "libmesh/meshfree_interpolation.h"

// C++ includes
#include <algorithm> // std::find
#include <cmath>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace libMesh
{
  // Forward Declarations
  template <typename T>
  class DenseVector;

  /**
   *
   */
  class PointInterpolationFunction : public FunctionBase<Number>
  {
  public:
    PointInterpolationFunction(const MeshfreeInterpolation &mfi,
                               Threads::spin_mutex &mutex);

    PointInterpolationFunction(const PointInterpolationFunction &);
    PointInterpolationFunction &operator=(const PointInterpolationFunction &);

    PointInterpolationFunction(PointInterpolationFunction &&) = default;
    PointInterpolationFunction &operator=(PointInterpolationFunction &&) = default;
    virtual ~PointInterpolationFunction() = default;

    virtual Number operator()(const Point &p,
                              const Real time = 0) override;

    /**
     * The actual initialization process.
     */
    void init() override;

    /**
     * Clears the function.
     */
    void clear() override;

    virtual std::unique_ptr<FunctionBase<Number>> clone() const override;

    virtual void operator()(const Point &p,
                            const Real time,
                            DenseVector<Number> &output) override;

  private:
    const MeshfreeInterpolation &_mfi;
    mutable std::vector<Point> _pts;
    mutable std::vector<Number> _vals;
    Threads::spin_mutex &_mutex;
  };

  /*----------------------- Inline functions ----------------------------------*/

  inline PointInterpolationFunction::PointInterpolationFunction(const MeshfreeInterpolation &mfi,
                                                                Threads::spin_mutex &mutex) : _mfi(mfi),
                                                                                              _mutex(mutex)
  {
  }

  inline PointInterpolationFunction &
  PointInterpolationFunction::operator=(const PointInterpolationFunction &other)
  {
    // Use copy-and-swap idiom
    PointInterpolationFunction tmp(other);
    std::swap(tmp, *this);
    return *this;
  }

  inline Number PointInterpolationFunction::operator()(const Point &p,
                                                       const Real /*time*/)
  {
    _pts.clear();
    _pts.push_back(p);
    _vals.resize(1);

    Threads::spin_mutex::scoped_lock lock(_mutex);

    _mfi.interpolate_field_data(_mfi.field_variables(), _pts, _vals);

    return _vals.front();
  }

  inline void PointInterpolationFunction::operator()(const Point &p,
                                                     const Real time,
                                                     DenseVector<Number> &output)
  {
    output.resize(1);
    output(0) = (*this)(p, time);
    return;
  }

  inline void PointInterpolationFunction::init()
  {
  }

  inline void PointInterpolationFunction::clear()
  {
  }

  inline std::unique_ptr<FunctionBase<Number>>
  PointInterpolationFunction::clone() const
  {
    return std::make_unique<PointInterpolationFunction>(_mfi, _mutex);
  }

} // namespace libMesh

#endif // LIBMESH_POINTS_INTERPOLATION_FUNCTION_H
