import glob
import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import convolve
from astropy.io import fits
import scipy.signal as sig
from numpy.linalg import lstsq
from decimal import Decimal
ccd_char='ccd_characteristics.txt'

def parseargs(argv=None):
    '''
    Options parser
    '''

    parser = argparse.ArgumentParser(description = "Data Analysis for Fe55")
    parser.add_argument("basename",nargs='?', type=str,action="store",
                              help="location of fits files")
    parser.add_argument("-o",dest = "ccd_char",action="store_true",default='ccd_characteristics.txt', help= "Output file")
    parser.add_argument("-wd",type=int, dest = "wd",default='50', help= "Window size of columns/rows to add estimate CTE")
    args = parser.parse_args(args=argv)

    if args.basename is None:
        msg = 'The base name was not provided'
        parser.error(msg)
    else:

        searchname = args.basename + '*/camra/*.fits'
        filenames = glob.glob(searchname)
        if not filenames:
            msg = 'No files found searching for: {:s}'.format(searchname)
            parser.error(msg)
        else:
            args.filenames = filenames
    return args


def main(argv=None):
    args = parseargs()
    ccd =['LL','LU','RL','RU']
    pdata,pcounts = {},{}
    gain,biastot = {},{}
    CTE = {'LL': [0,0], 'RL': [0,0],'LU':[0,0],'RU':[0,0]}
    ka_peak={'LL': [], 'RL': [],'LU':[],'RU':[]}
    pos={'LL': [], 'RL': [],'LU':[],'RU':[]}
    wd = args.wd

    for filename in args.filenames:
        position = filename[-11:-9]
        f = fits.open(filename)
        pos[position].append(f[0].data)
        f.close()
    for posi in pos:
        pos[posi] = np.array(pos[posi])

    for posicion in ccd:
        data,counts = np.unique(pos[posicion],return_counts=True)
        order=100
        peaks = sig.argrelextrema(np.log(np.convolve(counts,sig.gaussian(20, std=3),mode='same')),np.greater,order=order)[0]
        peaks = np.delete(peaks,np.argmax([counts[peak] for peak in peaks]))
        peak = peaks[np.argmax([counts[peak] for peak in peaks])]
        pdata[posicion],pcounts[posicion] = data[peak],counts[peak]
        plt.semilogy(pdata[posicion],pcounts[posicion],markersize=3,marker='o', color='r')
        plt.annotate('K_alpha', xy=(pdata[posicion], pcounts[posicion]))
        plt.title('CCD counts per DU')
        plt.ylabel('Counts')
        plt.xlabel('Signal,DU')
        plt.xlim((500,5000))
        plt.semilogy(data,counts)
        plt.savefig('plots/sig'+posicion+'.png')
        plt.show()
        print('Averaging over all pixels, K_alpha peak of %s in ND units should be: %4.1f'%(posicion,pdata[posicion]))
        print('Bias:%4.1f'%data[list(counts).index(np.amax(counts))])
        biastot[posicion] = data[list(counts).index(np.amax(counts))]
        #biastot = {'RL': 2328.7055833333334, 'RU': 2558.4873849999999, 'LL': 937.34572833333334, 'LU': 450.96781499999997}
        kalpha = pdata[posicion] - biastot[posicion]
        gain[posicion] = 1620./kalpha
        print('Peak minus offset is:%4.1f'%(kalpha))
        print('K_alpha photons produce 1620 electrons, the gain is therefore:%2.4f'%gain[posicion])
        print('RMS in electrons from the overscan is %3.3f'%(np.std(f[0].data[:,2070:2128])*gain[posicion]))
    
    for posicion in ccd:
        for row in xrange(len(pos['LL'][0,0,:2040])-wd):
            data,counts = np.unique(pos[posicion][:,:,row:row+wd],return_counts=True)
            bias = data[list(counts).index(np.amax(counts))]
            temp  = pdata[posicion]
            skip=False
            idx = (np.abs(data - temp)).argmin()
            if abs(temp - data[idx]) > 100 or len(counts)==0:
                continue
            data,counts = data[idx-100:idx+100],counts[idx-100:idx+100]
            counts = convolve(counts, sig.gaussian(wd,std=wd/3), mode='same')
            temppdata = data[np.argmax(counts)]
            ka_peak[posicion].append(temppdata)
    for posicion in ccd:
        plt.plot(ka_peak[posicion], markersize=1,marker='o', color='r', ls='')
        x=np.arange(len(ka_peak[posicion]))
        y=np.array(ka_peak[posicion])
        A = np.vstack([x, np.ones(len(x))]).T
        m, c = lstsq(A, y)[0]
        plt.plot(x, m*x + c, 'r', label='Fitted line')
        plt.title('K_alpha peaks, bias substracted per col')
        plt.ylabel('Signal, DU')
        plt.xlabel('column')
        plt.savefig('plots/col_sig'+posicion+'.png')
        plt.show()
        print len(ka_peak[posicion]),m,c
        CTE[posicion][0] = pow(1.+(len(ka_peak[posicion])*m)/c,1./len(ka_peak[posicion]))
        print('CTE per col %s in ND units is: %1.7f'%(posicion,CTE[posicion][0]))
    
    ka_peak={'LL': [], 'RL': [],'LU':[],'RU':[]}
    for posicion in ccd:
        for row in xrange(len(pos['LL'][0,:,0])-wd):
            data,counts = np.unique(pos[posicion][:,row:row+wd,:2040],return_counts=True)
            bias = data[list(counts).index(np.amax(counts))]
            temp  = pdata[posicion]
            skip=False
            idx = (np.abs(data - temp)).argmin()
            if abs(temp - data[idx]) > 100 or len(counts)==0:
                continue
            data,counts = data[idx-100:idx+100],counts[idx-100:idx+100]
            counts = np.convolve(counts,sig.gaussian(wd, std=wd/3),mode='same')
            temppdata = data[np.argmax(counts)]
            ka_peak[posicion].append(temppdata)
    for posicion in ccd:
        plt.plot(ka_peak[posicion], markersize = 1, marker = 'o', color = 'r', ls = '')
        x = np.arange(len(ka_peak[posicion]))
        y = np.array(ka_peak[posicion])
        A = np.vstack([x, np.ones(len(x))]).T
        m, c = lstsq(A, y)[0]
        plt.plot(x, m*x + c, 'r', label='Fitted line')
        plt.title('K_alpha peaks, bias substracted per row')
        plt.ylabel('Signal, DU')
        plt.xlabel('row')
        plt.savefig('plots/row_sig'+posicion+'.png')
        plt.show()
        print len(ka_peak[posicion]),m,c
        CTE[posicion][1] = pow(1.+(len(ka_peak[posicion])*m)/c,1./len(ka_peak[posicion]))
        print('CTE per row %s in ND units is: %1.7f'%(posicion,CTE[posicion][1]))

    output=open(args.ccd_char,'w')
    output.write('{0:14} {1:14} {2:14} {3:14} {4:14} \n'.format("#CCD","BIAS","GAIN","CTE(col)","CTE(row)"))
    for key in gain.keys():
        output.write('{0:14} {1:14} {2:14} {3:14} {4:14} \n'.format(key,biastot[key],gain[key],CTE[key][0],CTE[key][1]))
    output.close()
if __name__ == "__main__":
    main()
