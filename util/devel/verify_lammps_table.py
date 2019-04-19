#!/usr/bin/env python3

import sys

def read_table(filename):
    res = []
    with open(filename, 'r') as f:
        for line in f:
            items = line.split()
            if not len(line.strip()) or line[0] == '#' or len(items) != 4:
                continue
            res.append([float(x) for x in items[1:]])
    return res


def check_table(table):
    ferror = 0
    for i, values in enumerate(table[1:-1]):
        r = values[0]
        rprev = table[i][0]
        rnext = table[i+2][0]
        e = values[1]
        eprev = table[i][1]
        enext = table[i+2][1]
        f = values[2]
        fleft = -(e - eprev)/(r - rprev)
        fright = -(enext - e)/(rnext - r)
        if f < fleft and f < fright:
            print('FAIL @ {} < {} {} {}'.format(i + 1, fleft, f, fright))
            ferror += 1
        if f > fleft and f > fright:
            print('FAIL @ {} > {} {} {}'.format(i + 1, fleft, f, fright))
            ferror += 1
    print('Total errors: {}'.format(ferror))


def main(filename):
    if len(filename) != 1:
        print('Please provide a single filename as argument')
        sys.exit(1)
    check_table(read_table(filename[0]))


if __name__ == '__main__':
    main(sys.argv[1:])
