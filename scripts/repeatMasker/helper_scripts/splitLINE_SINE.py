#!/usr/bin/env python

import argparse 
#import re 


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='reformat the output of repeat masker for easy plotting in R')
    parser.add_argument('-i',type=str,dest='inpath',required=True,help="the last file i made ")
    parser.add_argument('-s',type=str,dest='sinepath',required=True,help="SINE-only tsv file being made")
    parser.add_argument('-l',type=str,dest='linepath',required=True,help="LINE-only tsv file being made")

    args = parser.parse_args()

#    inpath='/kyber/Data/Nanopore/Analysis/gilfunk/190821_brca1_circ_MOAR/watsHap/brca1_canu_assemblies/hap2/brca1_hap2.contigs.fasta.out'
#    outpath='/kyber/Data/Nanopore/Analysis/gilfunk/190821_brca1_circ_MOAR/watsHap/brca1_canu_assemblies/hap2_aluELEMENTS.tsv'

    infile=open( args.inpath  , 'r' )
    sinefile=open( args.sinepath , 'w')
    linefile=open( args.linepath , 'w')

    SINEs=[]
    LINEs=[]

    for line in infile: 
#        lmnts = []
        line=line.rstrip()
        tings=(line.split('\t'))
#        print('---------------')
#        for i in tings:
#            print( i ) 
        kind =  tings[4].split('/')
      #  print(kind[0])
        if kind[0] == 'SINE':
            SINEs.append(line)
        elif kind[0] == 'LINE':
            LINEs.append(line)
        else:
            pass
           
    print('TOTAL number of SINEs:')
    print(len(SINEs))
    print('TOTAL number of LINEs:')
    print(len(LINEs))

    for i in SINEs:
        sinefile.write(i+'\n')

    for i in LINEs:
        linefile.write(i+'\n')

    infile.close()
    sinefile.close()
    linefile.close()

#    print('this many types of alus')
#    print(len(kind_of_alus))
#print('hi')
