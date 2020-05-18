#!/usr/bin/env python

import argparse 
#import re 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='reformat the output of repeat masker for easy plotting in R')
    parser.add_argument('-i',type=str,dest='inpath',required=True,help="the output file from repeat mazkr")
    parser.add_argument('-o',type=str,dest='outpath',required=True,help="the new re-formatted tsv file being made")
    parser.add_argument('-a',type=str,dest='allele',required=True,help="which haplotype we wurkin on")
    args = parser.parse_args()

    infile=open( args.inpath  , 'r' )

    outfile=open( args.outpath , 'w')

    kind_of_alus=[]

    for line in infile: 
        lmnts = []
        line=line.rstrip()
        tings=(line.split(' '))
        for i in tings:
            if len(i) > 0:
                lmnts.append(i)
        if len(lmnts) > 0 :
            if lmnts[0] == 'SW':
                pass
            elif lmnts[0] == 'score':
                pass
            else:
                strand = '(+)'
                if lmnts[8] == 'C':
                    strand = '(-)'
      #      print( lmnts[5]+'\t'+lmnts[6]+'\t'+ lmnts[9]+'\t'+strand+'\n' )   
                kind_of_alus.append(lmnts[9])
                outfile.write(lmnts[4]+'\t'+ lmnts[5]+'\t'+lmnts[6]+'\t'+ lmnts[9]+'\t'+ lmnts[10] +'\t'+ strand+ '\t'+args.allele +   '\n' )

    infile.close()
    outfile.close()

#    print(len(kind_of_alus))
    kind_of_alus=set(kind_of_alus)

    print('this many types of alus')
    print(len(kind_of_alus))
#print('hi')
