#ifndef SCC_H_
#define SCC_H_

#include "averager.h"
// there is no need to have this header-only

class SCCCalculator {
public:
    SCCCalculator(const unsigned int max_k):
        _max_k(max_k), _n_stored(0)
    {
        _last_isis = new double[max_k - 1];
        memset(_last_isis, 0, (max_k - 1) * sizeof(double));
        _Iprod = new SimpleAverager[max_k];
    }

    virtual ~SCCCalculator()
    {
        delete[] _Iprod;
        delete[] _last_isis;
    }

    void feed_isi(const double isi)
    {
        _I.feed(isi);
        _Iprod[0].feed(isi * isi);
        for (int i = 0; i < _n_stored; i++)
            _Iprod[i + 1].feed(isi * _last_isis[i]);

        // isi storen, alle eins weiterruecken
        for (int i = _max_k - 2; i > 0; i--)
            _last_isis[i] = _last_isis[i - 1];
        _last_isis[0] = isi;
        if (_n_stored < _max_k - 1) _n_stored++;
    }

    void new_trial() {
         memset(_last_isis, 0, (_max_k - 1) * sizeof(double));
         _n_stored = 0;
    }

    double get_mean_I() const { return _I.get_mean(); }
    double get_var_I() const { return _I.get_var(); }
    double get_scc(const int k) const
    { 
        return (_Iprod[k].get_mean() - _I.get_mean() * _I.get_mean()) /
                _I.get_var();
    }
    double get_scc_sum() const
    {
        double res = 0;
        for (int i = 1; i < _max_k; i++)
        {
            res += get_scc(i);
        }
        return res;
    }


private:
    SimpleAverager* _Iprod;
    SimpleAverager _I;
    double* _last_isis;
    const int _max_k;
    int _n_stored;
        
};

#endif
