/* -*- c++ -*- */
/* 
 * Copyright 2017 <+YOU OR YOUR COMPANY+>.
 * 
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */


#ifndef INCLUDED_PHASMA_CLASSIFIER_H
#define INCLUDED_PHASMA_CLASSIFIER_H

#include <phasma/api.h>
#include <phasma/log.h>
#include <boost/circular_buffer.hpp>
#include <utility>

namespace gr {
  namespace phasma {

    /*!
     * \brief <+description+>
     *
     */
    class PHASMA_API classifier
    {
    public:
      classifier() {}

      ~classifier();

      void
      record_prediction(int predicted, int actual);

      float
      calc_pred_accuracy();

      void
      clear_pred_history();

      void
      calculate_confussion_matrix ();

      void
      print_confussion_matrix ();

    protected:
      classifier(size_t hist_size, size_t labels_num,
		 const std::vector<size_t> &labels);

      std::string
      decode_decision (int decision);

      boost::circular_buffer<std::pair <int, int>> d_pred_history;
      std::vector<std::vector<std::pair <int, float>>> d_confusion_matrix;

      std::vector<size_t> d_labels;

      size_t d_labels_num;
    };

  } // namespace phasma
} // namespace gr

#endif /* INCLUDED_PHASMA_CLASSIFIER_H */

