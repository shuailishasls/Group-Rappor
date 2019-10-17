#!/usr/bin/python2
# -*- coding: utf-8 -*

from __future__ import division
import sys
import optparse
import pandas as pd
import time
import rappor  # client library


def log(msg, *args):
    if args:
        msg = msg % args
    print(msg)


try:
    import fastrand
except ImportError:
    print("Native fastrand module not imported; see README for speedups")
    fastrand = None


def CreateOptionsParser():
    p = optparse.OptionParser()

    p.add_option(
        '--num-bits', type='int', metavar='INT', dest='num_bits', default=32,
        help='Number of bloom filter bits.')
    p.add_option(
        '--num-hashes', type='int', metavar='INT', dest='num_hashes', default=1,
        help='Number of hashes.')
    p.add_option(
        '--num-cohorts', type='int', metavar='INT', dest='num_cohorts', default=64,
        help='Number of cohorts.')

    p.add_option(
        '-p', type='float', metavar='FLOAT', dest='prob_p', default=0.25,
        help='Probability p')
    p.add_option(
        '-q', type='float', metavar='FLOAT', dest='prob_q', default=0.75,
        help='Probability q')
    p.add_option(
        '-f', type='float', metavar='FLOAT', dest='prob_f', default=0.5,
        help='Probability f')
    p.add_option(
        '--assoc-testdata', type='int', dest='assoc_testdata', default=0,
        help='Generate association testdata from true values on stdin.')

    choices = ['simple', 'fast']
    p.add_option(
        '-r', type='choice', metavar='STR',
        dest='random_mode', default='fast', choices=choices,
        help='Random algorithm (%s)' % '|'.join(choices))

    return p


def RapporClientSim(params, irr_rand, csv_in, row_count):
    """Read true values from csv_in and output encoded values on csv_out."""
    # header = ('client', 'cohort', 'bloom', 'prr', 'irr')
    # csv_out.writerow(header)

    # It would be more instructive/efficient to construct an encoder
    # instance up front per client, rather than one per row below.
    start_time = time.time()
    temp = []

    for i, (client_str, cohort_str, true_value) in enumerate(csv_in):
        try:
            cohort = int(cohort_str)
        except ValueError:
            pass

        secret = client_str
        enc = rappor.Encoder(params, cohort, secret, irr_rand)

        # Real users should call e.encode(). For testing purposes, we also want the PRR.
        bloom, prr, irr = enc._internal_encode(true_value)

        bloom_str = rappor.bit_string(bloom, params.num_bloombits)
        prr_str = rappor.bit_string(prr, params.num_bloombits)
        irr_str = rappor.bit_string(irr, params.num_bloombits)

        out_row = (client_str, cohort_str, '\t' + bloom_str, '\t' + prr_str, '\t' + irr_str)
        # csv_out.writerow(out_row)
        temp.append(out_row)

    return temp


def main(argv):
    (opts, argv) = CreateOptionsParser().parse_args(argv)

    # Copy flags into params
    try:
        with open('../csv/params.csv', 'rb') as csv_params:
            params = rappor.Params.from_csv(csv_params)
    except Exception as e:
        print(e)
        params = rappor.Params()
        params.num_bloombits = opts.num_bits
        params.num_hashes = opts.num_hashes
        params.num_cohorts = opts.num_cohorts
        params.prob_p = opts.prob_p
        params.prob_q = opts.prob_q
        params.prob_f = opts.prob_f

    if opts.random_mode == 'simple':
        irr_rand = rappor.SecureIrrRand(params)
    elif opts.random_mode == 'fast':
        if fastrand:
            log('Using fastrand extension')
            # NOTE: This doesn't take 'rand'.  It's seeded in C with srand().
            irr_rand = fastrand.FastIrrRand(params)
        else:
            log('Warning: fastrand module not importable; see README for build '
                'instructions.  Falling back to simple randomness.')
            irr_rand = rappor.SecureIrrRand(params)
    else:
        raise AssertionError
    # Other possible implementations:
    # - random.SystemRandom (probably uses /dev/urandom on Linux)
    # - HMAC-SHA256 with another secret?  This could match C++ byte for byte.
    #   - or srand(0) might do it.

    client_cohort = pd.read_csv('../csv/client_cohort_value.csv',
                                low_memory=False).values.tolist()
    result = RapporClientSim(params, irr_rand, client_cohort, len(client_cohort))
    pd.DataFrame(result).to_csv('../csv/reports.csv',
                                index=None, header=['client', 'cohort', 'bloom', 'prr', 'irr'])


if __name__ == "__main__":
    try:
        main(sys.argv)
    except RuntimeError as e:
        log('rappor_sim.py: FATAL: %s', e)
