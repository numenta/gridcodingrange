/* ---------------------------------------------------------------------
 * Numenta Platform for Intelligent Computing (NuPIC)
 * Copyright (C) 2019, Numenta, Inc.  Unless you have an agreement
 * with Numenta, Inc., for a separate license for this software code, the
 * following terms and conditions apply:
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero Public License version 3 as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Affero Public License for more details.
 *
 * You should have received a copy of the GNU Affero Public License
 * along with this program.  If not, see http://www.gnu.org/licenses.
 *
 * http://numenta.org/licenses/
 * ----------------------------------------------------------------------
 */

#ifndef NTA_BOX_EXPANSION
#define NTA_BOX_EXPANSION

#include <limits>
#include <vector>

/**
 * Expand a high-dimensional box in the positive direction to match a target
 * high-dimensional box. Both the starting and ending boxes start at the origin.
 * This class enumerates a set of box-shaped expansions that will expand the
 * first box into the second.
 */
class BoxExpansion {
public:
  template<typename It1, typename It2>
  void initialize(const It1 box_begin, const It1 box_end,
                  const It2 goal_begin, const It2 goal_end)
  {
    i_ = 0;
    progress_.assign(box_begin, box_end);
    goal_.assign(goal_begin, goal_end);
    ndim_ = std::distance(box_begin, box_end);
  }

  template<typename It>
  void setGoal(const It goal_begin, const It goal_end)
  {
    goal_.assign(goal_begin, goal_end);
    i_ = 0;
  }

  /**
   * Get the next box.
   *
   * Each box is an expansion along some dimension. This expansion box starts at
   * point [0, 0, ..., 0], except with one nonzero at dimension
   * nonzero_offset_dim.
   */
  bool getNext(size_t *nonzero_offset_dim,
               double *nonzero_offset_val,
               double shape[])
  {
    if (i_ >= ndim_) return false;

    const double growth = goal_[i_] - progress_[i_];
    if (growth > 0)
    {
      *nonzero_offset_dim = i_;
      *nonzero_offset_val = progress_[i_];
      std::copy(progress_.begin(), progress_.end(), shape);
      shape[i_] = growth;

      progress_[i_] = goal_[i_];
      i_++;
      return true;
    }
    else
    {
      // The larger box is not larger on this dimension. Skip to the next one.
      i_++;
      return this->getNext(nonzero_offset_dim, nonzero_offset_val, shape);
    }
  }

private:
  size_t i_;
  std::vector<double> progress_;
  std::vector<double> goal_;
  size_t ndim_;
};


class SelectiveIgnoranceBoxExpansion {
public:

  template<typename It>
  SelectiveIgnoranceBoxExpansion(
    const It scaledbox_begin, const It scaledbox_end,
    const It ignorebox_begin, const It ignorebox_end)
    : secondary_expanding_(false),
      scaledbox_(scaledbox_begin, scaledbox_end),
      ignorebox_(ignorebox_begin, ignorebox_end),
      ndim_(std::distance(scaledbox_begin, scaledbox_end))
  {
    baseline_factor_ = std::numeric_limits<double>::max();
    for (size_t i = 0; i < scaledbox_.size(); i++)
    {
      if (scaledbox_[i] > 0)
      {
        baseline_factor_ = std::min(baseline_factor_,
                                    ignorebox_[i] / scaledbox_[i]);
      }
    }

    expansion_factor_ = baseline_factor_ * 1.01;

    std::vector<double> initial(scaledbox_);
    std::vector<double> goal(scaledbox_);
    for (size_t i = 0; i < ndim_; i++)
    {
      initial[i] *= baseline_factor_;
      goal[i] *= expansion_factor_;
    }
    main_expansion_.initialize(
      initial.begin(), initial.end(),
      goal.begin(), goal.end());

    buffer_.assign(ndim_, -1);
  }

  void getNext(
    double offset[], double shape[], double *baseline_factor)
  {
    if (!secondary_expanding_)
    {
      if (main_expansion_.getNext(&main_nonzero_offset_dim_,
                                  &main_nonzero_offset_val_,
                                  shape))
      {
        if (main_nonzero_offset_val_ >=
            ignorebox_[main_nonzero_offset_dim_])
        {
          // The expansion region doesn't overlap with the ignore box. Simply
          // return the whole expansion region.
          std::fill(offset, offset + ndim_, 0);
          offset[main_nonzero_offset_dim_] =
            main_nonzero_offset_val_;
          *baseline_factor = baseline_factor_;
          return;
        }
        else
        {
          // The ignore-box overlaps with the expansion region. To check the
          // expansion region while ignoring the ignore-box, do a series of
          // checks that expand the ignore-box outward to fill the expansion
          // region.
          secondary_expanding_ = true;

          // The portion of the ignore box that falls within expansion region.
          // This is in the reference frame of the expansion region.
          buffer_ = ignorebox_;
          buffer_[main_nonzero_offset_dim_] =
            std::min(ignorebox_[main_nonzero_offset_dim_] -
                     main_nonzero_offset_val_, shape[main_nonzero_offset_dim_]);

          secondary_expansion_.initialize(buffer_.begin(), buffer_.end(),
                                          shape, shape + ndim_);
          return this->getNext(offset, shape, baseline_factor);
        }
      }
      else
      {
        // Expand and retry
        baseline_factor_ = expansion_factor_;
        expansion_factor_ *= 1.01;

        for (size_t i = 0; i < ndim_; i++)
        {
          buffer_[i] = scaledbox_[i]*expansion_factor_;
        }
        main_expansion_.setGoal(buffer_.begin(), buffer_.end());
        return this->getNext(offset, shape, baseline_factor);
      }
    }
    else
    {
      size_t nonzero_offset_dim;
      double nonzero_offset_val;
      if (secondary_expansion_.getNext(&nonzero_offset_dim,
                                       &nonzero_offset_val,
                                       shape))
      {
        std::fill(offset, offset + ndim_, 0);
        offset[nonzero_offset_dim] = nonzero_offset_val;

        // Convert to reference frame of the main expansion.
        offset[main_nonzero_offset_dim_] += main_nonzero_offset_val_;
        *baseline_factor = baseline_factor_;
        return;
      }
      else
      {
        secondary_expanding_ = false;
        return this->getNext(offset, shape, baseline_factor);
      }
    }
  }

private:
  BoxExpansion main_expansion_;
  size_t main_nonzero_offset_dim_;
  double main_nonzero_offset_val_;

  double baseline_factor_;
  double expansion_factor_;

  BoxExpansion secondary_expansion_;
  bool secondary_expanding_;

  std::vector<double> buffer_; // multi-use

  const std::vector<double> scaledbox_;
  const std::vector<double> ignorebox_;
  const size_t ndim_;
};


class MultiDirectionExpansion {
public:
  template<typename It>
  MultiDirectionExpansion(
    const It scaledbox_begin, const It scaledbox_end,
    const It ignorebox_begin, const It ignorebox_end,
    unsigned dimflags = (unsigned)-1)
    : bitvector_(0x0),
      started_(false),
      single_quadrant_expansion_(scaledbox_begin, scaledbox_end,
                                 ignorebox_begin, ignorebox_end),
      ndim_(std::distance(scaledbox_begin, scaledbox_end))
  {
    dimflags_ = (dimflags == (unsigned)-1)
      ? ((0x1 << ndim_) - 1)
      : dimflags;

    x0_unreflected_.assign(ndim_, -1);
    shape_.assign(ndim_, -1);
  }

  void getNext(double x0[], double shape[], double *baseline_factor)
  {
    if (!started_)
    {
      started_ = true;
    }
    else
    {
      while (true)
      {
        ++bitvector_;

        if (bitvector_ > dimflags_)
        {
          bitvector_ = 0x0;
          break;
        }

        if ((bitvector_ & ~dimflags_) == 0x0)
        {
          break;
        }
      }
    }

    if (bitvector_ == 0x0)
    {
      single_quadrant_expansion_.getNext(x0_unreflected_.data(),
                                         shape_.data(), &baseline_factor_);
    }

    // Perform appropriate reflection
    for (size_t i = 0; i < ndim_; i++)
    {
      const bool reflect = bitvector_ & (0x1 << i);
      if (reflect)
      {
        x0[i] = -(x0_unreflected_[i] + shape_[i]);
      }
      else
      {
        x0[i] = x0_unreflected_[i];
      }
    }
    std::copy(shape_.begin(), shape_.end(), shape);
    *baseline_factor = baseline_factor_;
  }

private:
  unsigned bitvector_;
  unsigned dimflags_;

  std::vector<double> x0_unreflected_;
  std::vector<double> shape_;
  double baseline_factor_;

  bool started_;

  SelectiveIgnoranceBoxExpansion single_quadrant_expansion_;

  const size_t ndim_;
};

#endif
