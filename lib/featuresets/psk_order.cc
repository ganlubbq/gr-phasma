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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <gnuradio/io_signature.h>
#include <phasma/featuresets/psk_order.h>
#include <volk/volk.h>

namespace gr
{
  namespace phasma
  {
    namespace featureset
    {

      psk_order::psk_order (size_t samples_num, float peak_threshold_db) :
	      d_samples_num (samples_num),
	      d_features_num (FEATURES_NUM),
	      d_peak_threshold_db (peak_threshold_db),
	      d_peak_threshold_lin (pow (10, peak_threshold_db / 10))
      {
	unsigned int alignment = volk_get_alignment ();

	d_outbuf = new float[FEATURES_NUM];
	d_inbuf = new gr_complex[d_samples_num];

	d_fftw_in_buf = (fftwf_complex *) fftwf_malloc (
	    d_samples_num * sizeof(fftw_complex));

	d_fftw_out_buf = (fftwf_complex *) fftwf_malloc (
	    d_samples_num * sizeof(fftwf_complex));

	d_shift = (gr_complex*) volk_malloc (d_samples_num * sizeof(gr_complex),
					     32);

	d_mag = (float*) volk_malloc (d_samples_num * sizeof(float), 32);

	// create fft plan to be used for channel power measurements
	d_fft = fftwf_plan_dft_1d (d_samples_num, d_fftw_in_buf, d_fftw_out_buf,
				  FFTW_FORWARD, FFTW_ESTIMATE);
      }

      psk_order::~psk_order ()
      {
	delete[] d_outbuf;
	delete[] d_inbuf;
	fftwf_destroy_plan (d_fft);
	fftwf_free (d_fftw_in_buf);
	fftwf_free (d_fftw_out_buf);
	volk_free (d_shift);
	volk_free (d_mag);
      }

      void
      psk_order::generate (const gr_complex* in)
      {
	memcpy (d_inbuf, in, d_samples_num * sizeof(gr_complex));
	d_outbuf[0] = multiply_and_detect ();

	memcpy (d_inbuf, d_fftw_in_buf, d_samples_num * sizeof(gr_complex));
	d_outbuf[1] = multiply_and_detect ();

	memcpy (d_inbuf, d_fftw_in_buf, d_samples_num * sizeof(gr_complex));
	d_outbuf[2] = multiply_and_detect ();
      }

      float
      psk_order::multiply_and_detect ()
      {
	volk_32fc_x2_multiply_32fc ((gr_complex*) d_fftw_in_buf, d_inbuf,
				    d_inbuf, d_samples_num);

	fftwf_execute (d_fft);

	/* Perform shifting and cropping on squared magnitude max noise-floor*/
	memcpy (&d_shift[0], &d_fftw_out_buf[d_samples_num / 2],
		sizeof(gr_complex) * (d_samples_num / 2));
	memcpy (&d_shift[d_samples_num / 2], d_fftw_out_buf,
		sizeof(gr_complex) * (d_samples_num / 2));

	volk_32fc_magnitude_32f (d_mag, (gr_complex*) d_fftw_out_buf,
				 d_samples_num);

	for (size_t i = 0; i < d_samples_num; i++) {
	  std::cout << 10 * log10 (d_mag[i]) << std::endl;
	}
	
	return 0;
      }

      size_t
      psk_order::get_features_num () const
      {
	return d_features_num;
      }

      float*
      psk_order::get_outbuf () const
      {
	return d_outbuf;
      }

      void
      psk_order::set_samples_num (size_t samples_num)
      {
	d_samples_num = samples_num;
      }

    } /* namespace featurest */
  } /* namespace phasma */
} /* namespace gr */

