// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeSurface.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Surfaces/ConeBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"
#include "Acts/Utilities/detail/RealQuadraticEquation.hpp"

namespace Acts {

/// @class ConeSurface
///
/// Class for a conical surface in the Tracking geometry.
/// It inherits from Surface.
///
/// The ConeSurface is special since no corresponding
/// Track parameters exist since they're numerical instable
/// at the tip of the cone.
/// Propagations to a cone surface will be returned in
/// curvilinear coordinates.

class ConeSurface : public Surface
{
public:
  ConeSurface() = delete;

  /// Constructor form HepTransform and an opening angle
  ///
  /// @param htrans is the transform to place to cone in a 3D frame
  /// @param alpha is the opening angle of the cone
  /// @param symmetric indicates if the cones are built to +/1 z
  ConeSurface(std::shared_ptr<const Transform3D> htrans,
              double                             alpha,
              bool                               symmetric = false);

  /// Constructor form HepTransform and an opening angle
  ///
  /// @param htrans is the transform that places the cone in the global frame
  /// @param alpha is the opening angle of the cone
  /// @param locZmin is the z range over which the cone spans
  /// @param locZmax is the z range over which the cone spans
  /// @param halfPhi is the openen angle for cone ssectors
  ConeSurface(std::shared_ptr<const Transform3D> htrans,
              double                             alpha,
              double                             locZmin,
              double                             locZmax,
              double                             halfPhi = M_PI);

  /// Constructor from HepTransform and ConeBounds
  ///
  /// @param htrans is the transform that places the cone in the global frame
  /// @param cbounds is the boundary class, the bounds must exit
  ConeSurface(std::shared_ptr<const Transform3D> htrans,
              std::shared_ptr<const ConeBounds>  cbounds);

  /// Copy constructor
  ///
  /// @param other is the source cone surface
  ConeSurface(const ConeSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param other is the source cone surface
  /// @param htrans is the additional transfrom applied after copying
  ConeSurface(const ConeSurface& other, const Transform3D& htrans);

  /// Constructor which accepts @c variant_data
  ///
  /// @param data the @c variant_data to build from
  ConeSurface(const variant_data& data);

  virtual ~ConeSurface();

  /// Assignment operator
  ///
  /// @param other is the source surface for the assignment
  ConeSurface&
  operator=(const ConeSurface& other);

  /// Implicit Constructor
  ///
  /// @param shift is the optional shift applied after cloning
  virtual ConeSurface*
  clone(const Transform3D* shift = nullptr) const final override;

  /// The binning position method - is overloaded for r-type binning
  ///
  /// @param bValue defines the type of binning applied in the global frame
  /// @return The return type is a vector for positioning in the global frame
  virtual const Vector3D
  binningPosition(BinningValue bValue) const final override;

  /// Return the surface type
  virtual SurfaceType
  type() const override;

  /// Return the measurement frame - this is needed for alignment, in particular
  ///  for StraightLine and Perigee Surface
  ///  - the default implementation is the the RotationMatrix3D of the transform
  ///
  /// @param gpos is the global position where the measurement frame is
  /// constructed
  /// @param mom is the momentum used for the measurement frame construction
  /// @return matrix that indicates the measurement frame
  const RotationMatrix3D
  referenceFrame(const Vector3D& gpos, const Vector3D& mom) const final;

  /// Return method for surface normal information
  ///
  /// @param lpos is the local position on the cone for which the normal vector
  /// is requested
  /// @return Vector3D normal vector in global frame
  const Vector3D
  normal(const Vector2D& lpos) const final;

  /// Return method for surface normal information
  ///
  /// @param gpos is the global position on the cone for which the normal vector
  /// is requested
  /// @return Vector3D normal vector in global frame
  const Vector3D
  normal(const Vector3D& gpos) const final;

  /// Normal vector return without argument
  using Surface::normal;

  // Return method for the rotational symmetry axis
  ///
  // @return This returns the local z axis
  virtual const Vector3D
  rotSymmetryAxis() const;

  /// This method returns the ConeBounds by reference
  virtual const ConeBounds&
  bounds() const final override;

  /// Local to global transformation
  ///
  /// @param lpos is the local position to be transformed
  /// @param mom is the global momentum (ignored in this operation)
  /// @param gpos is the global position shich is filled
  virtual void
  localToGlobal(const Vector2D& lpos,
                const Vector3D& mom,
                Vector3D&       gpos) const final override;

  /// Global to local transfomration
  ///
  /// @param gpos is the global position to be transformed
  /// @param mom is the global momentum (ignored in this operation)
  /// @param lpos is hte local position to be filled
  /// @return is a boolean indicating if the transformation succeeded
  virtual bool
  globalToLocal(const Vector3D& gpos,
                const Vector3D& mom,
                Vector2D&       lpos) const final override;

  /// @brief Straight line intersection schema - provides closest intersection
  /// and (signed) path length
  ///
  /// @param gpos The start position for the intersection
  /// @param gmom The start momentum for the intersection (will be normalized)
  /// @param navDir The navigation direction with respect to the momentum
  /// @param bcheck The boundary check to be used in this directive
  /// @param correct is an (optional) correction function pointer
  ///
  /// <b>mathematical motivation:</b>
  ///
  ///   The calculation will be done in the 3-dim frame of the cone,
  ///   i.e. the symmetry axis of the cone is the z-axis, x- and y-axis are
  /// perpendicular
  ///   to the the z-axis. In this frame the cone is centered around the origin.
  ///   Therefore the two points describing the line have to be first
  ///   recalculated
  /// into the new frame.
  ///   Suppose, this is done, the points of intersection can be
  ///   obtained as follows:<br>
  ///
  ///   The cone is described by the implicit equation
  ///   @f$x^2 + y^2 = z^2 \tan \alpha@f$
  ///   where @f$\alpha@f$ is opening half-angle of the cone  the and
  ///   the line by the parameter equation (with @f$t@f$ the
  ///   parameter and @f$x_1@f$ and @f$x_2@f$ are points on the line)
  ///   @f$(x,y,z) = \vec x_1 + (\vec x_2 - \vec x_2) t @f$.
  ///   The intersection is the given to the value of @f$t@f$ where
  ///   the @f$(x,y,z)@f$ coordinates of the line satisfy the implicit
  ///   equation of the cone. Inserting the expression for the points
  ///   on the line into the equation of the cone and rearranging to
  ///   the form of a  gives (letting @f$ \vec x_d = \frac{\vec x_2 - \vec
  ///   x_1}{|\vec x_2 - \vec x_1|} @f$):
  ///   @f$t^2 (x_d^2 + y_d^2 - z_d^2 \tan^2 \alpha) + 2 t (x_1 x_d +
  ///   y_1 y_d - z_1 z_d \tan^2 \alpha) + (x_1^2 + y_1^2 - z_1^2
  ///   \tan^2 \alpha) = 0 @f$
  ///   Solving the above for @f$t@f$ and putting the values into the
  ///   equation of the line gives the points of intersection. @f$t@f$
  ///   is also the length of the path, since we normalized @f$x_d@f$
  ///   to be unit length.
  ///
  /// @return is the Intersection object
  virtual Intersection
  intersectionEstimate(const Vector3D&      gpos,
                       const Vector3D&      gmom,
                       NavigationDirection  navDir,
                       const BoundaryCheck& bcheck = false,
                       CorrFnc correct = nullptr) const final override;

  /// the pathCorrection for derived classes with thickness
  ///
  /// @param gpos is the global potion at the correction point
  /// @param mom is the momentum at the correction point
  /// @return is the path correction due to incident angle
  virtual double
  pathCorrection(const Vector3D& gpos,
                 const Vector3D& mom) const final override;

  /// Return properly formatted class name for screen output
  virtual std::string
  name() const override;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  virtual variant_data
  toVariantData() const override;

protected:
  std::shared_ptr<const ConeBounds> m_bounds;  ///< bounds (shared)
};

#include "Acts/Surfaces/detail/ConeSurface.ipp"

}  // namespace