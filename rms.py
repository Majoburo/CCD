import glob
import argparse
import numpy as np
from astropy.io import fits






def parseargs(argv=None):
    '''
    Options parser
    '''

    parser = argparse.ArgumentParser(description = "Data Analysis for Fe55")
    parser.add_argument("basename",nargs='?', type=str,action="store",
                              help="location of fits files")
    args = parser.parse_args(args=argv)

    if args.basename is None:
        msg = 'The base name was not provided'
        parser.error(msg)
    else:

        args.fitsfiles = []
        searchname = args.basename + '*.fits'
        filenames = glob.glob(searchname)
        if not filenames:
            msg = 'No files found searching for: {:s}'.format(searchname)
            parser.error(msg)
        else:
            args.filenames = filenames
            args.fitsfiles = [fits.open(filenames[i])[0] for i in xrange(len(filenames))] 
    return args


def main(argv=None):
    args = parseargs()
    rms,mean={},{}
    for filename in args.filenames:
        position = filename[-11:-9]
        f = fits.open(filename)
        rms[position]=np.std(f[0].data[200:800,1000:2000])
        mean[position]=np.mean(f[0].data[200:800,1000:2000])
        f.close()
    print rms,mean

if __name__ == "__main__":
    main()
