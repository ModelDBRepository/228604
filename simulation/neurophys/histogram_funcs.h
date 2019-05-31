/* 
 * File:   increment_histogram.h
 * Author: fedro
 *
 * Created on May 27, 2015, 3:19 PM
 */

#ifndef HISTOGRAM_FUNCS_H
#define	HISTOGRAM_FUNCS_H

#include "dynamic_histogram.h"
#include "array.h"

namespace neurophys {

namespace hist {
    
// a little helper function to feed the increments of the time series x to a
// given histogram
    
void feed_array_increments_to_hist(const Array<double>& x, Histogram<unsigned int>* hist, 
        const double dt, const double tinc);
void feed_array_to_hist(const Array<double>& x, Histogram<unsigned int>* hist);

}
}

#endif	/* INCREMENT_HISTOGRAM_H */

