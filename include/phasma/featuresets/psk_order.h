/* -*- c++ -*- */
/* 
 * Copyright 2017 Kostis Triantafyllakis - ctriant.
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

#ifndef INCLUDED_PHASMA_PSK_ORDER_H
#define INCLUDED_PHASMA_PSK_ORDER_H

#include <phasma/api.h>
#include <fftw3.h>

#define FEATURES_NUM 3

namespace gr
{
  namespace phasma
  {
    namespace featureset
    {

    /*!
     * \brief <+description+>
     *
     */
    class PHASMA_API psk_order
    {
    public:
      psk_order(size_t samples_num, float peak_threshold_db);
      ~psk_order();

      void
      generate (const gr_complex* in);
      
      float
      multiply_and_detect ();

      size_t
      get_features_num () const;

      float*
      get_outbuf () const;

      void
      set_samples_num (size_t samples_num);

    private:
      float d_peak_threshold_db;
      float d_peak_threshold_lin;

      float* d_outbuf;
      float* d_mag;

      gr_complex* d_inbuf;
      gr_complex* d_shift;

      fftwf_complex* d_fftw_in_buf;
      fftwf_complex* d_fftw_out_buf;

      // fft plan for channel measurements
      fftwf_plan d_fft;

      size_t d_samples_num;
      size_t d_features_num;

    };
  }
} // namespace phasma
} // namespace gr

#endif /* INCLUDED_PHASMA_PSK_ORDER_H */

