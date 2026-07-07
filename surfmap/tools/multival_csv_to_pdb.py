#!/usr/bin/env python

import argparse, os, sys

# Script to convert multivalue output csv file to pdb file.


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required = True, help = "csv file from multivalue script (APBS suite)")
    parser.add_argument("-o", required = True, help = "output pdb file")
    args = parser.parse_args()

    ifile = args.i
    ofile = args.o

    indexline = 0

    outfile = open(ofile, "w+")

    with open(ifile) as f:
        for line in f:
            indexline += 1
            lsp = line.split(",")
            # numbers = [format(float(x), '.3f') for x in lsp]
            numbers = [float(x) for x in lsp]
            #print numbers
            value = numbers[3] if len(numbers) == 4 else 0.0

            # nb = '{0: >7}'.format(str(indexline))
            # x = '{0: >8}'.format(str(numbers[0]))
            # y = '{0: >8}'.format(str(numbers[1]))
            # z = '{0: >8}'.format(str(numbers[2]))
            #
            # if len(numbers) == 4:
            #     bf = '{0: >13}'.format(str(numbers[3]))
            # else:
            #     bf = '{0: >13}'.format("0.000")
            #
            # #print x,y,z,bf
            # stl = "ATOM" + nb + "  O   POS A   1    " + x + y + z + bf+ "\n"
            stl = (
                f"ATOM  {indexline:5d}  O   POS A   1    "
                f"{numbers[0]:8.3f}{numbers[1]:8.3f}{numbers[2]:8.3f}"
                f"{1.00:6.2f}{value:7.3f}          O\n"
            )
            #print stl
            outfile.write(stl)

    outfile.close()


if __name__ == "__main__":
    main()
